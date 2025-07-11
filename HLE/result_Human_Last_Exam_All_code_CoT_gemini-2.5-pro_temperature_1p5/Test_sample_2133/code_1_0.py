import numpy as np
from sympy import Matrix, permutedims
from math import sqrt
from itertools import permutations

def get_sn_character_table(n):
    """
    Returns the character table of the symmetric group S_n for n=4.
    The character table is represented as a dictionary mapping partition (character)
    to a dictionary mapping cycle structure (conjugacy class) to the character value.
    """
    if n == 4:
        # Conjugacy classes of S4, identified by cycle structure tuple
        # Partitions of 4 for rows (characters)
        # Partitions of 4 for columns (conjugacy classes)
        # Order of classes: (1,1,1,1), (2,1,1), (2,2), (3,1), (4)
        char_table = {
            "4":    {'1,1,1,1': 1, '2,1,1': 1, '2,2': 1, '3,1': 1, '4,1': 1},     # Trivial
            "3,1":  {'1,1,1,1': 3, '2,1,1': 1, '2,2': -1,'3,1': 0, '4,1': -1},    # Standard
            "2,2":  {'1,1,1,1': 2, '2,1,1': 0, '2,2': 2, '3,1': -1,'4,1': 0},
            "2,1,1":{'1,1,1,1': 3, '2,1,1': -1,'2,2': -1,'3,1': 0, '4,1': 1},     # Sign tensored with Standard
            "1,1,1,1":{'1,1,1,1': 1, '2,1,1': -1,'2,2': 1, '3,1': 1, '4,1': -1}   # Sign
        }
        return char_table
    else:
        raise NotImplementedError("Character table for S_n not available for n != 4")

def get_cycle_structure(perm):
    """Computes the cycle structure of a permutation."""
    n = len(perm)
    unvisited = list(range(n))
    cycles = []
    while unvisited:
        c = []
        i = unvisited[0]
        while i not in c:
            c.append(i)
            unvisited.remove(i)
            i = perm[i]
        cycles.append(len(c))
    return tuple(sorted(cycles, reverse=True))
    
def get_cycle_class_key(cycles):
    """Converts cycle structure to the key used in the character table."""
    key = ','.join(map(str, cycles))
    if key == '4': return '4,1' # Adjusting for my table format
    if key == '3': return '3,1'
    if key == '2': return '2,1,1'
    if key == '1': return '1,1,1,1'
    return key


def construct_mercer_matrix(n):
    """
    Constructs an n-nilpotent matrix with non-zero integer entries.
    """
    # Companion matrix for p(x) = x^n
    C = np.zeros((n, n), dtype=int)
    for i in range(1, n):
        C[i, i - 1] = 1

    # Vandermonde matrix for points 1, 2, ..., n
    S = np.array([[ (i+1)**j for j in range(n)] for i in range(n)], dtype=object)
    
    # Calculate S_inv and det(S) using sympy for precision
    S_sym = Matrix(S)
    S_inv_sym = S_sym.inv()
    det_S = S_sym.det()

    # M = det(S) * S * C * S^-1
    M_sym = det_S * S_sym * Matrix(C) * S_inv_sym
    
    # The result should have integer entries
    return np.array(M_sym.tolist()).astype(int)

def calculate_ratio(P):
    """
    Calculates the ratio of the logarithmic infinity norm to the Frobenius norm.
    """
    n = P.shape[0]
    # Frobenius norm
    frob_norm = np.linalg.norm(P, 'fro')
    if frob_norm == 0:
        return 0

    # Logarithmic infinity norm
    log_inf_norm = -np.inf
    for i in range(n):
        row_sum = P[i, i] + np.sum(np.abs(P[i, np.arange(n) != i]))
        if row_sum > log_inf_norm:
            log_inf_norm = row_sum
            
    return log_inf_norm / frob_norm

def calculate_largest_immanant(M):
    """
    Calculates the largest immanant of a matrix M.
    """
    n = M.shape[0]
    if n > 4:
        # Permutation-based calculation is too slow for n > 4.
        # This problem likely intends for a smaller n.
        return "Immanant calculation is too slow for n > 4."

    char_table = get_sn_character_table(n)
    
    # Map permutations to their cycle structure
    perms = list(permutations(range(n)))
    perm_prods = {p: np.prod([M[i, p[i]] for i in range(n)]) for p in perms}
    
    # Group permutations by conjugacy class
    class_sums = {class_key: 0 for class_key in char_table['4'].keys()}
    for p in perms:
        cycles = get_cycle_structure(p)
        key = get_cycle_class_key(cycles)
        class_sums[key] += perm_prods[p]

    # Calculate immanants
    immanants = {}
    for char_name, char_values in char_table.items():
        imm = sum(char_values[class_key] * class_sums[class_key] for class_key in class_sums)
        immanants[char_name] = imm
        
    largest_immanant = max(immanants.values())
    return largest_immanant

def solve_for_n(n):
    """
    Main function to execute the plan for a given n.
    """
    print(f"Solving for n = {n}")
    
    # 1. Construct Mercer Matrix
    M_n = construct_mercer_matrix(n)
    print("Constructed Mercer Matrix M_n:")
    print(M_n)

    # 2. Compute Popov Normal Form (RREF)
    P_n_sym, _ = Matrix(M_n).rref()
    # The Popov form may require pivots to be positive, RREF from sympy is sufficient here.
    P_n = np.array(P_n_sym.tolist()).astype(float)
    print("\nPopov-like Form (RREF) P_n:")
    print(P_n)

    # 3. Calculate Ratio
    ratio = calculate_ratio(P_n)
    print(f"\nRatio of Logarithmic Infinity Norm to Frobenius Norm: {ratio}")

    # 4. Find the Largest Immanant
    largest_imm = calculate_largest_immanant(M_n)
    print(f"Largest Immanant of M_{n}: {largest_imm}")
    return largest_imm

if __name__ == '__main__':
    # Set n to a value for which the calculation is feasible, e.g., n=4.
    n_val = 4
    final_answer = solve_for_n(n_val)
    print(f"\nFinal answer for the largest immanant for n={n_val}.")
    # The final answer format is specified in the problem.
    # The code calculates this value and prints it at the end.
    print(f"<<<{final_answer}>>>")