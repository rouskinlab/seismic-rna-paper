
lines = list()
with open("out-de-lajarte-2024/menrich.csv") as f:
    for line in f:
        if line.count(",") != 5:
            raise ValueError(line)
        lines.append(line[line.index(",") + 1:])

text = "".join(lines)
print(text)

