import collections

def solve():
    """
    Calculates the minimum number of initially infected sites for an n x n grid.
    The problem reduces to finding the minimum number of sites (1s) to place
    in a grid such that no 2x2 subgrid of healthy sites (0s) remains.

    A pattern that achieves this is to infect sites (i,j) where (i + 2j) % 4 == C.
    We calculate the number of infected sites for each C in {0,1,2,3} and find the minimum.
    """
    n = 14
    m = 4

    # n_k stores the number of indices i in {0, ..., n-1} such that i % m == k
    n_k = collections.defaultdict(int)
    for i in range(n):
        n_k[i % m] += 1
    
    print(f"For n={n} and modulus m={m}, the counts of indices for each residue class are:")
    for k in sorted(n_k.keys()):
        print(f"Number of indices i where i % {m} == {k}: {n_k[k]}")
    print("-" * 20)

    min_sites = n * n
    
    # We will demonstrate the calculation for C=0
    # The equation is: count = sum_{j_mod_m} (n_{j_mod_m} * n_{(C - 2*j_mod_m)_mod_m})
    
    counts = []
    for C in range(m):
        current_count = 0
        eq_parts = []
        # Calculate count for each C
        # i = C - 2j (mod m)
        for j_mod_m in range(m):
            # Number of j's where j % m == j_mod_m
            num_j = n_k[j_mod_m]
            # Required residue class for i
            i_mod_m = (C - 2 * j_mod_m) % m
            # Number of i's with that residue class
            num_i = n_k[i_mod_m]
            term = num_j * num_i
            current_count += term
            eq_parts.append(f"{num_j}*{num_i}")
        
        counts.append(current_count)
        print(f"For C = {C}:")
        print(f"Infected sites = {' + '.join(eq_parts)} = {current_count}")
        
    min_sites = min(counts)

    print("-" * 20)
    print(f"The minimum number of sites to infect for n={n} is: {min_sites}")

solve()
<<<49>>>