import os
import subprocess
import sys
from itertools import combinations, product
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from seismicrna.core.arg import opt_mask_polya
from seismicrna.core.logs import logger, set_config
from seismicrna.core.random import stochastic_round
from seismicrna.core.rel import MATCH, RelPattern
from seismicrna.core.rna import from_ct
from seismicrna.core.seq import parse_fasta
from seismicrna.core.task import as_list_of_tuples, dispatch
from seismicrna.sim import (fastq as fastq_mod,
                            fold as fold_mod,
                            params as params_mod,
                            ref as ref_mod)
from seismicrna.sim.muts import load_pmut, calc_pmut_pattern

rng = np.random.default_rng()


# The longest amplicon possible with 150x150 nt paired-end reads and
# clipping 4 nt from each end (during relate) is 292 nt.
AMPLICON_MAX_LENGTH = 292  # nt
PRIMER_LENGTH = 20  # nt
MEAN_MUT_PER_READ = 2.0
UNPAIRED_PAIRED_MUT_RATIO = 3.0
# Simulations in test_avg_paired.py show that the average fraction of
# paired bases that RNAstructure predicts is about 58%.
MEAN_F_PAIRED = 0.58
SIM_DIR = Path("sim")
SAMPLES_DIR = SIM_DIR.joinpath("samples")


def sim_params(reflen: int,
               k: int,
               refnum: int,
               num_cpus: int,
               max_ident: float = 2/3,
               max_r: float = 0.5**0.5,
               max_tries: int = 200):
    assert reflen >= 1
    assert k >= 1
    assert refnum >= 0
    ref = f"rna-{reflen}-{k}-{refnum}"
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
            if reflen <= AMPLICON_MAX_LENGTH:
                # Use amplicon mode.
                read_length_mean = reflen
                length_fmean = read_length_mean / reflen
                length_fvar = 0.
                primer_length = PRIMER_LENGTH
            else:
                # Use fragment mode.
                read_length_mean = 250  # nt
                length_fmean = read_length_mean / reflen
                read_length_stdev = 10  # nt
                length_fvar = ((read_length_stdev / reflen) ** 2
                               / (length_fmean * (1 - length_fmean)))
                # Fragment mode has no primers.
                primer_length = 0
            # Make the average number of mutations per read 2, and make
            # the unpaired mutation rate 3 times the paired mutation rate,
            # which is about the ratio seen in experimental data.
            # First, estimate the number of bases in the read that have
            # the ability to mutate, which is the read length minus both
            # primers (if any) times 0.5 (because on average 50% of the
            # bases are As and Cs).
            num_mutable_bases = (read_length_mean - 2 * primer_length) * 0.5
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
                           length_fmean=length_fmean,
                           length_fvar=length_fvar,
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
            # Give every position in a primer binding site 0 mutations.
            pmut_file = ct_file.with_suffix(".muts.csv")
            primer_positions = np.concatenate(
                    [list(range(primer_length)),
                     list(range(reflen - primer_length, reflen))]
            ) + 1
            if primer_positions.size > 0:
                pmut = pd.read_csv(pmut_file,
                                   index_col=[0, 1],
                                   header=[0, 1, 2])
                pmut.loc[primer_positions] = 0.
                pmut.loc[primer_positions, str(MATCH)] = 1.
                pmut.to_csv(pmut_file)
            # For cluster calculations, mask out the bases that will not
            # be used.
            region = structures[0].region
            region.mask_gu()
            region.mask_polya(opt_mask_polya.default)
            region.mask_list(primer_positions)
            # Ensure every cluster is sufficiently different from every
            # other cluster, as measured by the Pearson correlation over
            # A and C bases.
            mus = calc_pmut_pattern(load_pmut(pmut_file),
                                    RelPattern.from_counts(count_del=True))
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


SIMS = list(range(60))
REFLENS = [140, 280, 420]
KS = [1, 2, 3, 4]
NUMS_READS = [5000,
              10000,
              20000,
              50000,
              100000,
              200000]
NUM_READS_POP = {
    0: 200000,
    2: 15000,
    3: 60000,
    4: 100000,
    5: 50000,
    6: 40000,
    7: 30000,
    8: 20000
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
    "nobias": {0: 1.0},
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
                     reflen: int,
                     k: int,
                     sim: int):
    proportions = PROP_READS_SAMP[sample_type]
    proportions_array = np.array(list(proportions.values()))
    assert np.isclose(proportions_array.sum(), 1.)
    nums_reads = stochastic_round(num_reads * proportions_array,
                                  preserve_sum=True)
    sample = f"sample-{sample_type}-{num_reads}"
    ref = f"rna-{reflen}-{k}-{sim}"
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


def jackpot(n: int, r: int):
    sample = np.arange(n)
    for _ in range(r):
        sample = rng.choice(sample, size=n, replace=True)
    return sample


def sim_fastq_jackpot(fastq_prefix: Path, r: int):
    print("Jackpot", fastq_prefix)
    fastq_prefix_name = fastq_prefix.name
    sample_in = fastq_prefix.parent.name
    if sample_in.count("-biased-") != 1:
        raise ValueError(sample_in)
    sample_out = sample_in.replace("-biased-", f"-jackpot-{r}-")
    sample = None
    for mate in [1, 2]:
        fastq_name = f"{fastq_prefix_name}_R{mate}.fq"
        fastq_file_in = fastq_prefix.with_name(fastq_name)
        fastq_file_out = fastq_prefix.parent.parent.joinpath(sample_out, fastq_name)
        if fastq_file_out.is_file():
            continue
        reads = list()
        with open(fastq_file_in) as f:
            for line in f:
                read = [line]
                read.extend(f.readline() for _ in range(3))
                reads.append("".join(read))
        if sample is None:
            sample = jackpot(len(reads), r)
            fastq_file_out.parent.mkdir(exist_ok=True)
        with open(fastq_file_out, "x") as f:
            for i, x in enumerate(sample):
                # Append a unique ID number to each read name because
                # some reads will otherwise have duplicate names.
                f.write(reads[x].replace("\n", f"-{i}\n", 1))


def simulate(num_cpus: int):
    # Simulate reference sequences and structures.
    fastas_cts = dispatch(sim_params,
                          num_cpus=num_cpus,
                          pass_num_cpus=True,
                          raise_on_error=False,
                          as_list=True,
                          ordered=True,
                          args=product(REFLENS, KS, SIMS))
    fastas = list()
    cts = list()
    for fasta, ct in fastas_cts:
        fastas.append(fasta)
        cts.append(ct)
    # Merge all reference sequences into one FASTA file.
    combined_fasta = Path("all-rnas.fa").resolve()
    with open(combined_fasta, "w") as fc:
        for indiv_fasta in fastas:
            with open(indiv_fasta) as fi:
                fc.write(fi.read())
    # Write a regions file.
    regions_file = Path("regions.csv").resolve()
    with open(regions_file, "w") as f:
        f.write("Region,Reference,5' End,3' End,Forward Primer,Reverse Primer\n")
        for ref in parse_fasta(combined_fasta, None):
            _, reflen_str, k_str, refnum_str = ref.split("-")
            reflen = int(reflen_str)
            if reflen <= AMPLICON_MAX_LENGTH:
                f.write(",".join(map(str, [f"amp-{reflen}",
                                           ref,
                                           PRIMER_LENGTH + 1,
                                           reflen - PRIMER_LENGTH,
                                           "",
                                           "\n"])))
    return
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
                          REFLENS,
                          KS,
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
    # Simulate the jackpotted FASTQ files.
    RS = [0, 2, 4, 6, 8, 10, 12]
    dispatch(sim_fastq_jackpot,
             num_cpus=num_cpus,
             pass_num_cpus=False,
             raise_on_error=True,
             as_list=True,
             ordered=False,
             args=[(SAMPLES_DIR.joinpath("sample-biased-200000",
                                         f"rna-280-1-{sim}"),
                    r)
                   for sim, r in product(SIMS, RS)])


if __name__ == "__main__":
    set_config(verbosity=4,
               exit_on_error=True)
    num_cpus = int(sys.argv[1])
    simulate(num_cpus)
