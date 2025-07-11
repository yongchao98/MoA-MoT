import math
from itertools import permutations
import numpy as np

def get_permutation_cycle_type(p):
    """
    Finds the cycle type of a permutation.
    The permutation p is a tuple, e.g., (0, 1, 3, 2).
    Returns a tuple of cycle lengths, sorted in descending order.
    """
    n = len(p)
    visited = [False] * n
    cycles = []
    for i in range(n):
        if not visited[i]:
            cycle = []
            j = i
            while not visited[j]:
                visited[j] = True
                cycle.append(j)
                j = p[j]
            cycles.append(len(cycle))
    return tuple(sorted(cycles, reverse=True))

def get_s4_char_table_and_maps():
    """
    Provides the character table for S_4 and returns a list of dictionaries,
    each mapping a cycle type to its character value for one irreducible representation.
    """
    # Conjugacy classes of S4, represented by partitions of 4: (1^4), (2,1^2), (2,2), (3,1), (4)
    partitions = [(1, 1, 1, 1), (2, 1, 1), (2, 2), (3, 1), (4,)]
    
    # Character table for S4 from theory
    # Rows correspond to irreps: X[4], X[3,1], X[2,2], X[2,1,1], X[1^4]
    table = [
        [1, 1, 1, 1, 1],       # Trivial (associated with Permanent)
        [3, 1, -1, 0, -1],     # Standard representation
        [2, 0, 2, -1, 0],
        [3, -1, -1, 0, 1],
        [1, -1, 1, 1, -1]      # Sign (associated with Determinant)
    ]

    char_maps = []
    for i in range(len(table)):
        char_map = {partitions[j]: table[i][j] for j in range(len(partitions))}
        char_maps.append(char_map)
        
    return char_maps

def build_Mn(n):
    """
    Constructs the specific n x n n-nilpotent matrix with non-zero integer entries.
    This construction yields a matrix with row sums equal to zero.
    """
    if n < 2:
        return np.array([])
    M = np.zeros((n, n), dtype=int)
    
    # Set the first row
    M[0, :n-1] = 1
    M[0, n-1] = -(n - 1)
    
    # Recursively define subsequent rows up to n-1
    for i in range(1, n - 1):
        M[i] = M[i-1].copy()
        M[i, i] += 1
        M[i, n-1] -= 1
        
    # The last row is identical to the (n-1)-th row
    M[n-1] = M[n-2].copy()
    
    return M

def calculate_immanant(matrix, char_map):
    """
    Calculates the immanant of a matrix for a given character map.
    """
    n = matrix.shape[0]
    total = 0
    
    # Iterate through all permutations in S_n
    for p in permutations(range(n)):
        
        # Calculate the product term for the permutation
        prod = 1
        for i in range(n):
            prod *= matrix[i, p[i]]
            
        # Get character value for the permutation's cycle type
        cycle_type = get_permutation_cycle_type(p)
        chi_val = char_map[cycle_type]
        
        total += chi_val * prod
        
    return total

def solve_for_n(n):
    """
    Solves the problem for a given n.
    """
    print(f"For n={n}, the specific Mercer matrix M_n is:")
    Mn = build_Mn(n)
    print(Mn)

    # Note: Character table and immanant calculation gets very slow for n > 8
    # The following part is specific to n=4.
    if n != 4:
        print(f"\nImmanant calculation is implemented only for n=4.")
        return

    char_maps = get_s4_char_table_and_maps()
    immanants = []
    
    char_names = ['[4] (Permanent)', '[3,1]', '[2,2]', '[2,1,1]', '[1,1,1,1] (Determinant)']
    
    print("\nCalculating immanants for M_4:")
    for i, char_map in enumerate(char_maps):
        imm = calculate_immanant(Mn, char_map)
        immanants.append(imm)
        print(f"imm_{char_names[i]}:\t{imm}")

    largest_immanant = max(immanants)
    
    print(f"\nThe largest immanant for M_4 is {largest_immanant}.")
    return largest_immanant

# --- Main Execution ---
# We solve for n=4 as a representative example.
solve_for_n(4)
