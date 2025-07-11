import collections

def get_schubert_cell_counts(k, n):
    """
    Calculates the number of Schubert cells N_p for each dimension p.
    Cells are indexed by partitions (l_1, ..., l_k) with n-k >= l_1 >= ... >= l_k >= 0.
    The dimension of the cell is p = sum(l_i).
    """
    m = n - k
    dim = k * m
    counts = collections.defaultdict(int)
    
    # Recursive function to find partitions
    def find_partitions(p_idx, current_sum, max_val, partition):
        if p_idx == k:
            counts[current_sum] += 1
            return
            
        for val in range(max_val, -1, -1):
            find_partitions(p_idx + 1, current_sum + val, val, partition + [val])

    find_partitions(0, 0, m, [])
    
    # Convert defaultdict to a list for dimensions 0 to dim
    N = [counts[p] for p in range(dim + 1)]
    return N

def solve_torsion_rank():
    """
    Calculates the rank of the torsion subgroup of the integral cohomology
    ring of G(3,5).
    """
    k = 3
    n = 5
    dim = k * (n - k)

    # Step 1: Get Schubert cell counts N_p
    N = get_schubert_cell_counts(k, n)
    print(f"The number of Schubert cells N_p for p=0..{dim} are: {N}")

    # Step 2: Define Betti numbers b_p
    # For G(3,5), b_0=1, b_4=1, others are 0.
    b = [0] * (dim + 1)
    b[0] = 1
    b[4] = 1
    print(f"The Betti numbers b_p are: {b}")

    # Step 3: Recursively calculate torsion ranks t_p for homology
    t = [0] * (dim + 1)
    t_prev = 0
    print("\nCalculating the homology torsion ranks t_p using the formula t_p = N_p - b_p - t_{p-1}:")
    for p in range(dim + 1):
        t[p] = N[p] - b[p] - t_prev
        print(f"p={p}: t_{p} = N_{p} - b_{p} - t_{p-1} = {N[p]} - {b[p]} - {t_prev} = {t[p]}")
        t_prev = t[p]

    print(f"\nThe ranks of the homology torsion groups t_p are: {t}")
    
    # Step 4: Calculate the total rank of the cohomology torsion subgroup
    # Rank of Tors(H^p) is t_{p-1}. Total rank is sum over p>=1.
    # This is equivalent to sum of t_p for p>=0.
    total_torsion_rank = sum(t)
    
    print("\nThe rank of the torsion subgroup of the integral cohomology ring is the sum of the ranks of the homology torsion groups.")
    
    equation_parts = [str(val) for val in t]
    print(f"Total Rank = {' + '.join(equation_parts)} = {total_torsion_rank}")
    
    return total_torsion_rank

if __name__ == '__main__':
    rank = solve_torsion_rank()
    # The final answer is directly printed by the function.
