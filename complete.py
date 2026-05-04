import numpy as np
from scipy.special import factorial
from itertools import permutations

def constrained_permutations_fast(n):
    """
    Build permutations directly without filtering.
    For each pair (i, i+1), treat them as a unit and only emit
    orderings where i precedes i+1.
    
    Equivalent to: take all permutations, then for each pair,
    fix the relative order — giving n! / 2^(n//2) results.
    """
    
    pairs = list(range(0, n - 1, 2))   # start of each constrained pair
    lone  = [n - 1] if n % 2 != 0 else []  # unpaired last element if n is odd

    # Permute over "slots" and assign pair order deterministically
    elements = list(range(n))

    def _build(chosen, remaining):
        if not remaining:
            yield tuple(chosen)
            return
        seen = set()
        for i, val in enumerate(remaining):
            if val in seen:
                continue
            seen.add(val)
            # If val is the larger of a pair (odd, and val-1 is even), skip —
            # it must come after its partner, which hasn't been placed yet
            if val % 2 == 1 and (val - 1) not in chosen:
                continue
            rest = remaining[:i] + remaining[i+1:]
            yield from _build(chosen + [val], rest)

    yield from _build([], elements)

def count_crossings(values):
    """
    values: flat list [a0, b0, a1, b1, ...] where (values[2i], values[2i+1]) are pairs.
    Counts pairs of intervals that properly cross, in O(n log n)
    where n is the number of pairs.
    """
    # Build intervals and sort by left endpoint
    intervals = [(values[i], values[i+1]) for i in range(0, len(values), 2)]
    intervals.sort(key=lambda x: (x[0], -x[1]))

    # Count inversions in the right endpoints via merge sort
    rights = [b for _, b in intervals]
    _, count = _merge_sort(rights)
    return count

def _merge_sort(arr):
    if len(arr) <= 1:
        return arr, 0
    mid = len(arr) // 2
    left,  lc = _merge_sort(arr[:mid])
    right, rc = _merge_sort(arr[mid:])
    merged, sc = _merge(left, right)
    return merged, lc + rc + sc

def _merge(left, right):
    result, inversions, i, j = [], 0, 0, 0
    while i < len(left) and j < len(right):
        if left[i] <= right[j]:
            result.append(left[i]); i += 1
        else:
            inversions += len(left) - i
            result.append(right[j]); j += 1
    result.extend(left[i:])
    result.extend(right[j:])
    return result, inversions

def count_crossings_alt(values):
    result = 0
    for i, val in enumerate(values):
        result += max(0, i-val)
    return result

max_graph_size = 10
graph_sizes = range(2, max_graph_size + 1)
for size in graph_sizes:
    tot_num = factorial(2*size)/(2**size)
    tot_good = 0
    for perm in constrained_permutations_fast(2*size):
        tot_good += 1 if count_crossings_alt(perm)%2 == 0 else 0
    print(tot_good, tot_num)
    print(2*(tot_good-tot_num/2)/tot_num)


# 0 2 3 1 4 5


