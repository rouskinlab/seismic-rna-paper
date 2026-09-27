import os
import subprocess
import sys
from collections import defaultdict
from itertools import combinations, product
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from seismicrna.core.arg import opt_mask_polya
from seismicrna.core.logs import logger, set_config
from seismicrna.core.random import stochastic_round
from seismicrna.core.rel import MATCH, RelPattern
from seismicrna.core.rna import from_ct, db_to_ct
from seismicrna.core.seq import DNA, parse_fasta
from seismicrna.core.task import as_list_of_tuples, dispatch
from seismicrna.sim import (fastq as fastq_mod,
                            fold as fold_mod,
                            params as params_mod,
                            ref as ref_mod)
from seismicrna.sim.muts import load_pmut, calc_pmut_pattern

rng = np.random.default_rng()


MEAN_MUT_PER_READ = 3.0
UNPAIRED_PAIRED_MUT_RATIO = 3.0
# Simulations in test_avg_paired.py show that the average fraction of
# paired bases that RNAstructure predicts is about 58%.
MEAN_F_PAIRED = 0.58
SIM_DIR = Path("sim")
SAMPLES_DIR = SIM_DIR.joinpath("samples")
READ_LENGTH_MEAN = 240  # nt


def sim_params(reflen: int,
               k: int,
               refnum: int,
               region: int,
               num_cpus: int,
               max_ident: float = 1/2,
               max_r: float = 0.5,
               max_tries: int = 200):
    assert reflen >= 1
    assert k >= 1
    assert refnum >= 0
    ref = f"long-rna-{refnum}-region-{region}"
    for i in range(max_tries):
        force = bool(i)
        # Simulate the reference sequence.
        fasta = ref_mod.run(refs=ref, ref=ref, reflen=reflen, force=force)
        # Simulate the structures.
        ct_file, = fold_mod.run(fasta,
                                fold_max=k,
                                num_cpus=num_cpus,
                                force=force)
        # Ensure the number of simulated structures is correct.
        structures = list(from_ct(ct_file))
        if len(structures) == k:
            # Make the average number of mutations per read 2, and make
            # the unpaired mutation rate 3 times the paired mutation rate,
            # which is about the ratio seen in experimental data.
            # First, estimate the number of bases in the read that have
            # the ability to mutate, which is the read length minus both
            # primers (if any) times 0.5 (because on average 50% of the
            # bases are As and Cs).
            num_mutable_bases = READ_LENGTH_MEAN * 0.5
            # Then estimate the number of paired and unpaired bases.
            num_paired_bases = num_mutable_bases * MEAN_F_PAIRED
            num_unpaired_bases = num_mutable_bases - num_paired_bases
            # The average number of mutations per read must equal the
            # number of mutable paired bases times the paired mutation
            # rate plus the number of mutable unpaired bases times the
            # unpaired mutation rate:
            # MEAN_MUT_PER_READ = (num_paired_bases * m_paired
            #                      +
            #                      num_unpaired_bases * m_unpaired)
            #                   = (num_paired_bases * m_paired
            #                      +
            #                      (num_unpaired_bases
            #                       * (UNPAIRED_PAIRED_MUT_RATIO * m_paired)))
            # m_paired = (MEAN_MUT_PER_READ /
            #             (num_paired_bases + (num_unpaired_bases
            #                                  * UNPAIRED_PAIRED_MUT_RATIO))
            m_paired = (MEAN_MUT_PER_READ /
                        (num_paired_bases + (num_unpaired_bases
                                             * UNPAIRED_PAIRED_MUT_RATIO)))
            m_unpaired = UNPAIRED_PAIRED_MUT_RATIO * m_paired
            params_mod.run(ct_file=[ct_file],
                           # All clusters have the same proportion.
                           clust_conc=1.e18,
                           # These parameters come from the output of
                           # parameterize/de-lajarte-2024.sh,
                           # except the mutation rates (am and cm).
                           vmut_paired=0.009146551012491061,
                           vmut_unpaired=0.01664274411629216,
                           pmut_paired=[('loq', 0.002157742662795042),
                                        ('am', m_paired),
                                        ('ac', 0.16728751343904066),
                                        ('ag', 0.3414862605146254),
                                        ('at', 0.42750935133585743),
                                        ('cm', m_paired),
                                        ('ca', 0.2647246848457093),
                                        ('cg', 0.1893355357454597),
                                        ('ct', 0.5289675777148524)],
                           pmut_unpaired=[('loq', 0.00353998550800204),
                                          ('am', m_unpaired),
                                          ('ac', 0.15950268981311605),
                                          ('ag', 0.2699218482634623),
                                          ('at', 0.5305362802833817),
                                          ('cm', m_unpaired),
                                          ('ca', 0.3039153717793367),
                                          ('cg', 0.16144348919771792),
                                          ('ct', 0.5171874972937138)],
                           num_cpus=num_cpus,
                           force=force)
            # For cluster calculations, mask out the bases that will not
            # be used.
            pmut_file = ct_file.with_suffix(".muts.csv")
            region = structures[0].region
            region.mask_gu()
            region.mask_polya(opt_mask_polya.default)
            # Ensure every cluster is sufficiently different from every
            # other cluster, as measured by the Pearson correlation over
            # A and C bases.
            mus = calc_pmut_pattern(load_pmut(pmut_file),
                                    RelPattern.from_counts(count_del=False))
            mus = mus.loc[region.unmasked_int]
            sufficiently_diff = True
            for col1, col2 in combinations(mus.columns, 2):
                mus1 = mus.loc[:, col1]
                mus2 = mus.loc[:, col2]
                fraction_ident = np.mean(np.isclose(mus1, mus2))
                if fraction_ident >= max_ident:
                    sufficiently_diff = False
                pearson_r = np.corrcoef(mus1, mus2)[0, 1]
                if pearson_r >= max_r:
                    sufficiently_diff = False
            if sufficiently_diff:
                mus_file = pmut_file.with_name("mus.csv")
                mus.to_csv(mus_file)
                logger.action(f"Wrote mutation rates to {mus_file}")
                return fasta, ct_file
    raise RuntimeError(f"Failed to simulate reference {ref}")


def concat_params(refnum: int, regions: int):
    refseqs = list()
    db_strings = list()
    pmuts = list()
    for region in range(regions):
        ref = f"long-rna-{refnum}-region-{region}"
        # Read the reference sequence.
        fasta = SIM_DIR.joinpath("refs", f"{ref}.fa")
        fasta_refs = list(parse_fasta(fasta, DNA))
        if len(fasta_refs) != 1:
            raise ValueError(f"Expected 1 reference in {fasta}, got {len(fasta_refs)}")
        fasta_ref, fasta_refseq = fasta_refs[0]
        if fasta_ref != ref:
            raise ValueError(f"Expected reference {ref} in {fasta}, got {fasta_ref}")
        refseqs.append(fasta_refseq)
        # Read the structure.
        params_dir = SIM_DIR.joinpath("params", ref, "full")
        ct_file = params_dir.joinpath("simulated.ct")
        structures = list(from_ct(ct_file))
        k = len(structures)
        db_strings.append([structure.db_structure for structure in structures])
        # Read the mutation rates.
        pmut_file = params_dir.joinpath("simulated.muts.csv")
        pmut = pd.read_csv(pmut_file, index_col=[0, 1], header=[0, 1, 2])
        start = 1 + sum(map(len, refseqs[:-1]))
        pmut.index = pd.MultiIndex.from_arrays([np.arange(start, start + pmut.index.size),
                                                pmut.index.get_level_values("Base")],
                                               names=pmut.index.names)
        pmuts.append(pmut)
    # Generate the combined reference sequence and structure.
    ref = f"long-rna-{refnum}"
    fasta = SIM_DIR.joinpath("refs", f"{ref}.fa")
    # Concatenate the reference sequences.
    refseq = DNA("".join(map(str, refseqs)))
    reflen = len(refseq)
    with open(fasta, "w") as f:
        f.write(f">{ref}\n{refseq}")
    params_dir = SIM_DIR.joinpath("params", ref, "full")
    params_dir.mkdir(exist_ok=True, parents=True)
    # List all possible structures.
    structs_nums = list(product(*[range(len(db_string)) for db_string in db_strings]))
    pmut = pd.DataFrame(
        np.nan,
        index=pd.MultiIndex.from_arrays([np.arange(1, reflen + 1),
                                        list(refseq)],
                                        names=["Position", "Base"]),
        columns=pd.MultiIndex.from_product(
            [[2,16,32,64,113,128,177,209,225,241,1],
                [len(structs_nums)],
                np.arange(1, len(structs_nums) + 1)],
            names=["Relationship", "K", "Cluster"],
        )
    )
    db_file = params_dir.joinpath(f"simulated.db")
    with open(db_file, "w") as f:
        pass
    for cluster, struct_nums in enumerate(structs_nums, start=1):
        db_string = "".join(db_strings[region][struct_num]
                            for region, struct_num in enumerate(struct_nums))
        with open(db_file, "a") as f:
            f.write(f">structure-{cluster}\n")
            if cluster == 1:
                f.write(f"{refseq.tr()}\n")
            f.write(f"{db_string}\n")
        for region, struct_num in enumerate(struct_nums):
            logger.detail(f"Assigning pmut for region {region}, struct_num {struct_num}, cluster {cluster}")
            pmut.loc[pmuts[region].index, pmut.columns.get_level_values("Cluster") == cluster] = pmuts[region].loc[:, pmuts[region].columns.get_level_values("Cluster") == str(struct_num + 1)].values
    pmut_file = params_dir.joinpath("simulated.muts.csv")
    pmut.to_csv(pmut_file)
    mus = calc_pmut_pattern(pmut, RelPattern.from_counts(count_del=False))
    mus.to_csv(params_dir.joinpath("mus.csv"))
    ct_file = db_to_ct(db_file)
    length_fmean = READ_LENGTH_MEAN / reflen
    read_length_stdev = 10  # nt
    length_fvar = ((read_length_stdev / reflen) ** 2
                    / (length_fmean * (1 - length_fmean)))
    params_mod.run(ct_file=[ct_file],
                    length_fmean=length_fmean,
                    length_fvar=length_fvar,
                    # All clusters have the same proportion.
                    clust_conc=1.e18)
    return ct_file


def sim_fastq_pop(ct_file: Path, min_mut_gap: int, num_reads: int, num_cpus: int):
    return fastq_mod.run(input_path=[],
                         param_dir=[ct_file.parent],
                         sample=f"population-{min_mut_gap}",
                         paired_end=True,
                         read_length=150,
                         min_mut_gap=min_mut_gap,
                         num_reads=num_reads,
                         fq_gzip=False,
                         num_cpus=num_cpus)


SIMS = list(range(10))
NUMS_READS = [1000000]
NUM_READS_POP = {
    0: 181527,
    2:  60000,
    3: 250000,
    4: 450000,
    5: 180000,
    6: 120000,
    7: 100000,
    8:  30000,
}
# The proportion of each population going into the biased sample was
# calculated based on the data in parameterize/out-morandi-2021/pooled/
# graph_mutdist/sars-cov-2/full/mutdist_filtered_m.csv
#
# Gap Mutated Mutated-NULL Mut/Mut-NULL Normalized Difference
#   0  664225  1798502.607     0.200175   0.181527   0.181527
#   1  267855  1338102.838     0.200175   0.181527   0.0
#   2  291084  1147433.897     0.253682   0.230049   0.048523
#   3  441972   920700.664     0.480038   0.435318   0.205269
#   4  726853   861399.434     0.843804   0.765196   0.329878
#   5  846775   886808.246     0.954856   0.865903   0.100707
#   6  881504   850348.352     1.036638   0.940066   0.074163
#   7  888662   816303.479     1.088641   0.987224   0.047158
#   8  950564   862010.306     1.102729   1.0        0.012775
PROP_READS_SAMP = {
    "biased": {0: 0.181527,
               1: 0.000000,
               2: 0.048523,
               3: 0.205269,
               4: 0.329878,
               5: 0.100707,
               6: 0.074163,
               7: 0.047158,
               8: 0.012775}
}
LINES_PER_READ = 4


def count_fastq_reads(fastq: str | Path):
    num_lines = int(subprocess.check_output(["wc", "-l", fastq]).split()[0])
    num_reads, rem = divmod(num_lines, LINES_PER_READ)
    if rem:
        raise ValueError(f"{fastq} has {num_lines} lines")
    return num_reads


def iter_fastq_reads(fastq: str | Path, read_nums: Iterable[int]):
    with open(fastq) as f:
        position = 0
        for read_num in read_nums:
            if position > read_num:
                raise ValueError("read_nums is not sorted or contains repeats "
                                 "or values that are not non-negative integers")
            while position < read_num:
                # Skip lines in groups of 4.
                for _ in range(LINES_PER_READ):
                    f.readline()
                position += 1
            for _ in range(LINES_PER_READ):
                line = f.readline()
                if not line:
                    raise ValueError("No more lines")
                yield line
            position += 1


def sim_fastq_sample(sample_type: str,
                     num_reads: int,
                     sim: int):
    proportions = PROP_READS_SAMP[sample_type]
    proportions_array = np.array(list(proportions.values()))
    assert np.isclose(proportions_array.sum(), 1.)
    nums_reads = stochastic_round(num_reads * proportions_array,
                                  preserve_sum=True)
    sample = f"sample-{sample_type}-{num_reads}"
    ref = f"long-rna-{sim}"
    read_nums = dict()
    for mate in [1, 2]:
        fastq_name = f"{ref}_R{mate}.fq"
        fastq_sample = SAMPLES_DIR.joinpath(sample, fastq_name)
        logger.routine(f"Began writing {fastq_sample}")
        fastq_sample.parent.mkdir(exist_ok=True, parents=False)
        with open(fastq_sample, "w") as f:
            for pop, num_reads_from_pop in zip(proportions,
                                               nums_reads,
                                               strict=True):
                if num_reads_from_pop <= 0:
                    continue
                population = f"population-{pop}"
                fastq_population = SAMPLES_DIR.joinpath(population, fastq_name)
                pop_read_nums = read_nums.get(pop)
                if pop_read_nums is None:
                    pop_num_reads = count_fastq_reads(fastq_population)
                    logger.detail(f"Choosing {num_reads_from_pop} reads from {fastq_population} "
                                  f"(containing {pop_num_reads} reads) into {fastq_sample}")
                    pop_read_nums = np.sort(rng.choice(pop_num_reads,
                                                       num_reads_from_pop,
                                                       replace=False,
                                                       shuffle=False))
                    read_nums[pop] = pop_read_nums
                for line in iter_fastq_reads(fastq_population, pop_read_nums):
                    if line.startswith("@"):
                        line = line.replace("@", f"@{population}-")
                    f.write(line)
        logger.routine(f"Ended writing {fastq_sample}")
    return SAMPLES_DIR.joinpath(sample, ref)


REGIONS = [(200, 1),
           (400, 2),
           ( 50, 1),
           (150, 3),
           (100, 2),
           (100, 2),
           (200, 1)]
SIMS = list(range(60))


def simulate(num_cpus: int):
    # Simulate reference sequences and structures.
    fastas_cts = dispatch(sim_params,
                          num_cpus=num_cpus,
                          pass_num_cpus=True,
                          raise_on_error=False,
                          as_list=True,
                          ordered=True,
                          args=[(reflen, k, sim, region)
                                for region, (reflen, k) in enumerate(REGIONS)
                                for sim in SIMS])
    cts = list()
    for sim in SIMS:
        cts.append(concat_params(sim, len(REGIONS)))
    # Simulate FASTQ files of populations of reads with varying amounts
    # of bias.
    fastqs_pop = list()
    for min_mut_gap, num_reads in NUM_READS_POP.items():
        for fastq_pair in dispatch(sim_fastq_pop,
                                   num_cpus=num_cpus,
                                   pass_num_cpus=True,
                                   args=as_list_of_tuples(cts),
                                   raise_on_error=True,
                                   as_list=False,
                                   ordered=True,
                                   kwargs=dict(min_mut_gap=min_mut_gap,
                                               num_reads=num_reads)):
            fastqs_pop.extend(fastq_pair)
    # Simulate FASTQ files of samples from the populations.
    dispatch(sim_fastq_sample,
             num_cpus=num_cpus,
             pass_num_cpus=False,
             raise_on_error=True,
             as_list=True,
             ordered=False,
             args=product(PROP_READS_SAMP,
                          NUMS_READS,
                          SIMS))
    # Delete the population FASTQ files.
    fastqs_pop_dirs = set()
    for fastq in fastqs_pop:
        fastqs_pop_dirs.add(fastq.parent)
        fastq.unlink()
        logger.action(f"Deleted {fastq}")
    for fastqs_pop_dir in fastqs_pop_dirs:
        fastqs_pop_dir.rmdir()
        logger.action(f"Deleted {fastqs_pop_dir}")


if __name__ == "__main__":
    set_config(verbosity=4,
               exit_on_error=True)
    num_cpus = int(sys.argv[1])
    simulate(num_cpus)
