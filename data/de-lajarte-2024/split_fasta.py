# Split the FASTA file of all references into one FASTA file per plate.

import os 

from seismicrna.core.seq import parse_fasta, DNA

PREFIX = "mRNA_plate-"

# Numbers of the references (0-indexed) in each 384-well plate.
PLATE_TO_REFS = {
    1: list(range(0, 384)),
    2: list(range(384, 768)),
    3: list(range(768, 1152)),
    # Plate 4 is split into two chunks.
    4: list(range(3066, 3258)) + list(range(4406, 4551)),
    # Plate 5 has 378 references, not 384.
    5: list(range(1536, 1914)),
    # At this point, the indexes get shifted so that each plate starts
    # at 6 less than a multiple of 384, rather than an even multiple.
    6: list(range(1914, 2298)),
    7: list(range(2298, 2682)),
    8: list(range(2682, 3066)),
    9: list(range(3258, 3642)),
    # Plates 10 and 11 have 382 references, not 384.
    10: list(range(3642, 4024)),
    11: list(range(4024, 4406)),
    # Plate 12 is in the position where plate 4 should be.
    12: list(range(1152, 1536)),
}

REF_TO_PLATE = dict()
for platenum, refnums in PLATE_TO_REFS.items():
    for refnum in refnums:
        if refnum in REF_TO_PLATE:
            raise ValueError(refnum)
        REF_TO_PLATE[refnum] = platenum


def get_plate_fasta(platenum: int, prefix: str):
    return f"{prefix}{platenum}.fa"


def split_fasta(fasta: str, prefix: str):
    # Erase the FASTA file for each plate.
    for platenum in PLATE_TO_REFS:
        with open(get_plate_fasta(platenum, prefix), "w") as f:
            pass
    for refnum, (refname, refseq) in enumerate(parse_fasta(fasta, DNA)):
        platenum = REF_TO_PLATE[refnum]
        with open(get_plate_fasta(platenum, prefix), "a") as f:
            f.write(f">{refname}\n{refseq}\n")
        print(f"Wrote {refname} to {platenum}")


if __name__ == "__main__":
    split_fasta("amplicon_rnas.fa", PREFIX)
    #fastas = [get_plate_fasta(platenum, PREFIX) for platenum in PLATE_TO_REFS]
    #os.system(f"cat {' '.join(fastas)} > mRNA_plates.fa")

