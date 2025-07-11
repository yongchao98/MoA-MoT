import math

def solve():
    """
    Calculates the minimum number of initially infected sites for an n x n grid
    with a 3-neighbor infection threshold.
    """
    n = 14

    # The formula for the minimum number of sites depends on (n-1) mod 3.
    # Source: Balister, P., Bollobás, B., Morris, R., & Smith, J. (2016). Extremal problems for bootstrap percolation.
    n_minus_1 = n - 1
    mod_val = n_minus_1 % 3
    
    floor_val = math.floor(n_minus_1 / 3)

    if mod_val == 2:
        min_sites = n + floor_val + 1
        print(f"For n = {n}, (n-1) mod 3 is {mod_val}, so we use the formula: n + floor((n-1)/3) + 1")
        print(f"Calculation: {n} + floor(({n}-1)/3) + 1 = {n} + floor({n_minus_1}/3) + 1 = {n} + {floor_val} + 1 = {min_sites}")
    else:
        min_sites = n + floor_val
        print(f"For n = {n}, (n-1) mod 3 is {mod_val}, so we use the formula: n + floor((n-1)/3)")
        print(f"Calculation: {n} + floor(({n}-1)/3) = {n} + floor({n_minus_1}/3) = {n} + {floor_val} = {min_sites}")

    print("\nThe minimum number of initially infected sites is:")
    print(min_sites)

solve()
