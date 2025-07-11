# Python code to calculate the result
from functools import lru_cache

@lru_cache(maxsize=None)
def solve(target, dims):
    if not dims:
        return 1 if target == 0 else 0
    
    d, m, rest_dims = dims[0], dims[1], dims[2:]
    count = 0
    for i in range(target // d + 1):
        # Number of ways to choose multiplicities for this dimension is i+m-1 choose m-1
        # For m=2, this is i+1. For m=1, this is 1.
        if m == 2:
            term_count = i + 1
        else: # m=1
            term_count = 1
        
        count += term_count * solve(target - i * d, rest_dims)
    return count

# dims_mults = [(6,1), (5,2), (4,2), (1,2)]
# result = solve(1000, tuple(x for y in dims_mults for x in y))
# The above recursive approach is slow. A dynamic programming approach is better.

memo = {}
def count_partitions(n, parts_mults):
    if (n, parts_mults) in memo:
        return memo[(n, parts_mults)]
    if n == 0:
        return 1
    if n < 0 or not parts_mults:
        return 0
    
    d, m = parts_mults[0]
    rest_parts = parts_mults[1:]
    
    res = 0
    for i in range(n // d + 1):
        # Number of ways to partition i*d into m parts of size d
        # is the number of ways to partition i into m parts of size 1
        # which is number of solutions to x_1+...+x_m = i, which is C(i+m-1, m-1)
        from math import comb
        term_count = comb(i + m - 1, m - 1)
        res += term_count * count_partitions(n - i * d, rest_parts)
        
    memo[(n, parts_mults)] = res
    return res

# parts_with_multiplicities = ((6, 1), (5, 2), (4, 2), (1, 2))
# result = count_partitions(1000, parts_with_multiplicities)
# This calculation yields the final answer.