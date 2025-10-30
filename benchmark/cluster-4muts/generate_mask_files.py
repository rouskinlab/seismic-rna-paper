import sys


REVCOMP = str.maketrans({"A": "T", "C": "G", "G": "C", "T": "A"})

def revcomp(seq: str):
    return seq[::-1].translate(REVCOMP)


def read_primer_lengths(regions_file: str):
    lengths = dict()
    with open(regions_file) as f:
        f.readline()  # Skip the header.
        for line in f:
            amp, ref, end5, end3, _ = line.split(",", maxsplit=4)
            reflen = int(ref.split("-")[1])
            fwdlen = int(end5) - 1
            revlen = reflen - int(end3)
            lengths[ref] = (fwdlen, revlen)
    return lengths


def generate_rf_mask_line(ref: str, fwd: str, rev: str):
    masks = list()
    if fwd:
        masks.append(fwd)
    if rev:
        masks.append(rev)
    return f"{','.join([ref] + masks)}\n" if masks else ""


def generate_sm_primers_line(ref: str, fwd: str, rev: str):
    return f">{ref}\n{fwd} {revcomp(rev)}\n" if (fwd and rev) else ""


def generate_files(rf_mask_file: str,
                   sm_primers_file: str,
                   fasta_file: str,
                   regions_file: str):
    lengths = read_primer_lengths(regions_file)
    rf_lines = list()
    sm_lines = list()
    with open(fasta_file) as f:
        for line in f:
            if not line.startswith(">"):
                raise ValueError(line)
            ref = line.rstrip()[1:]
            seq = f.readline().rstrip()
            fwdlen, revlen = lengths.get(ref, (0, 0))
            fwdseq = seq[:fwdlen] if fwdlen > 0 else ""
            revseq = seq[-revlen:] if revlen > 0 else ""
            rf_lines.append(generate_rf_mask_line(ref, fwdseq, revseq))
            sm_lines.append(generate_sm_primers_line(ref, fwdseq, revseq))
    with open(rf_mask_file, "w") as f:
        f.write("".join(rf_lines))
    with open(sm_primers_file, "w") as f:
        f.write("".join(sm_lines))


if __name__ == "__main__":
    generate_files(*sys.argv[1:])

