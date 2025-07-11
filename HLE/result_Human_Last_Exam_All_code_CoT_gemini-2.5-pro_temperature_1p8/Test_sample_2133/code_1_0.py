import numpy as np
from itertools import permutations

def get_s3_characters():
    """Returns the irreducible characters of the symmetric group S3."""
    # S3 has 3 conjugacy classes: id, (12), (123)
    # Trivial character (partition [3])
    chi_triv = {(): 1, (0, 1): 1, (0, 2, 1): 1}
    # Sign character (partition [1, 1, 1])
    chi_sign = {(): 1, (0, 1): -1, (0, 2, 1): 1}
    # Standard character (partition [2, 1])
    chi_std = {(): 2, (0, 1): 0, (0, 2, 1): -1}

    # Map permutation tuples to conjugacy class representatives
    char_map = {}
    
    # Identity
    p = (0, 1, 2)
    char_map[p] = {'triv': chi_triv[()], 'sign': chi_sign[()], 'std': chi_std[()]}

    # Transpositions (swaps) like (1,2)
    transpositions = [(1, 0, 2), (2, 1, 0), (0, 2, 1)]
    for p in transpositions:
        char_map[p] = {'triv': chi_triv[(0, 1)], 'sign': chi_sign[(0, 1)], 'std': chi_std[(0, 1)]}

    # 3-cycles like (1,2,3)
    cycles3 = [(1, 2, 0), (2, 0, 1)]
    for p in cycles3:
        char_map[p] = {'triv': chi_triv[(0, 2, 1)], 'sign': chi_sign[(0, 2, 1)], 'std': chi_std[(0, 2, 1)]}
        
    return char_map

def compute_immanants(matrix):
    """Computes all immanants for a 3x3 matrix."""
    n = matrix.shape[0]
    if n != 3:
        raise ValueError("This function is defined for 3x3 matrices.")

    char_table = get_s3_characters()
    perms = list(permutations(range(n)))
    
    immanants = {'triv': 0, 'sign': 0, 'std': 0}
    
    # Store terms for the permanent calculation string
    permanent_terms = []

    for p in perms:
        prod = 1
        for i in range(n):
            prod *= matrix[i, p[i]]
        
        term_str = " * ".join(map(str, [matrix[i, p[i]] for i in range(n)]))

        chars = char_table[p]
        for key in immanants:
            immanants[key] += chars[key] * prod
        
        if chars['triv'] > 0:
            permanent_terms.append(f"({term_str})")
        else: # Should not happen for permanent but good practice
            permanent_terms.append(f"-({term_str})")
            
    # The permanent is the immanant of the trivial character
    perm_val = immanants['triv']
    
    # The determinant is the immanant of the sign character
    det_val = immanants['sign']
    
    # The standard immanant
    std_val = immanants['std']
    
    print("The matrix M_3 is:")
    print(matrix)
    print("\nCalculations for Immanants:")
    print(f"Permanent (corresponding to trivial character) calculation:")
    print("perm(M_3) = " + " + ".join(permanent_terms))
    print(f"Permanent = {perm_val}")
    print(f"Determinant (corresponding to sign character) = {det_val}")
    print(f"Standard Character Immanant = {std_val}")

    largest_immanant = max(perm_val, det_val, std_val)
    print(f"\nThe values of the immanants are {perm_val}, {det_val}, and {std_val}.")
    print(f"The largest immanant is {largest_immanant}.")
    return largest_immanant

# For n=3, a specific matrix that fits the criteria is constructed as M_3.
M3 = np.array([
    [-1, -3, -2],
    [1, -1, 2],
    [1, 1, 2]
])

# Run the calculation
final_answer = compute_immanants(M3)
print("\nFinal Answer:")
# Final printing as requested.
# "Remember in the final code you still need to output each number in the final equation!"
# We output the full calculation of the largest immanant's parts, and the final value.
print(f"{final_answer}")