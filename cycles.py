from functools import lru_cache
from fractions import Fraction
from math import factorial

def cycle_neighbors(n):
    nbrs = [set() for _ in range(n)]
    if n >= 2:
        for i in range(n):
            j = (i + 1) % n
            nbrs[i].add(j)
            nbrs[j].add(i)
    return nbrs

def exact_b_cycle(n):
    nbrs = cycle_neighbors(n)
    total_words = factorial(2 * n) // (2 ** n)   # number of double-occurrence words

    @lru_cache(None)
    def dp(started_mask, open_order):
        open_order = list(open_order)

        # all chords have been opened and closed
        if started_mask == (1 << n) - 1 and not open_order:
            return 1

        total = 0

        # option 1: open a new chord
        for v in range(n):
            if not (started_mask >> v) & 1:
                total += dp(started_mask | (1 << v), tuple(open_order + [v]))

        # option 2: close an existing open chord
        for i, v in enumerate(open_order):
            later_open = open_order[i + 1:]

            # when v closes, it crosses exactly the neighbors that are still open
            # and were opened after v
            crossings = sum(1 for u in later_open if u in nbrs[v])
            sign = -1 if crossings % 2 else 1

            new_open = tuple(open_order[:i] + open_order[i + 1:])
            total += sign * dp(started_mask, new_open)

        return total

    numerator = dp(0, ())
    return Fraction(numerator, total_words)

for n in range(1, 8):
    value = exact_b_cycle(n)
    print(f"b(C_{n}) = {value} = {float(value):.12f}")
