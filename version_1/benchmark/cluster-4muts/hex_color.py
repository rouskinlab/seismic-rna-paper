import sys


def decimal_to_hex(n: int):
    if not isinstance(n, int):
        raise TypeError(n)
    if not 0 <= n < 256:
        raise ValueError(n)
    return f"{hex(n)[2:]:0>2}"


def triplet_to_hex_color(r: int, g: int, b: int):
    return f"#{''.join(map(decimal_to_hex, [r, g, b]))}"


if __name__ == "__main__":
    r, g, b = map(int, sys.argv[1:])
    print(triplet_to_hex_color(r, g, b))
    
