import numpy as np
from sympy.combinatorics.permutations import Permutation
from sympy.combinatorics.representation import Representation

def calculate_largest_immanant(M):
    """
    Calculates the largest immanant for a given square matrix M.
    """
    n = M.shape[0]
    
    # S_n is the symmetric group of degree n
    # In sympy, it's PermutationGroup(SymmetricGroup(n))
    # However, we only need the characters, which can be generated directly.
    # For S_3, the characters are well-known and can be hardcoded
    # or generated. We will generate them for generality.
    
    # Get all irreducible representations for S_n
    # For S_3, these correspond to partitions [[3]], [[2, 1]], [[1, 1, 1]]
    reps = Representation.irreps(n)
    
    immanants = []
    
    print(f"The matrix M_n (for n=3) is:\n{M}\n")

    # A dictionary to hold characters for S_n
    # chi[permutation] = character_value
    characters = {}
    for r in reps:
        # r is a representation object
        # r.name gives the partition
        # r(p) gives the matrix for permutation p
        # character is the trace of this matrix
        partition_name = r.name
        char_table = {p: r(p).trace() for p in r.group}
        characters[partition_name] = char_table

    print("Calculating immanants for all irreducible characters of S_3:")
    
    # Get all permutations of S_n
    perms = list(Permutation.group_elements(n))

    for partition, char_table in characters.items():
        imm = 0
        for p in perms:
            # p is a Permutation object
            # p.array_form gives the one-line notation [p(0), p(1), ...]
            # We need to convert to 1-based for matrix indexing
            perm_list = [i + 1 for i in p.array_form]
            
            prod = 1
            for i in range(n):
                # M[i, j] where j = p(i)
                prod *= M[i, perm_list[i]-1]
            
            # Get character value for permutation p
            chi_val = char_table[p]
            imm += chi_val * prod
        
        immanants.append(imm)
        print(f"  Immanant for partition {partition}: {imm}")

    largest_immanant = max(immanants)
    return largest_immanant

# For n=3, the specific matrix M_3 that satisfies the conditions
M3 = np.array([
    [1, 1, 2],
    [1, -3, 2],
    [1, -5, 2]
])

# Calculate and print the largest immanant
largest_imm = calculate_largest_immanant(M3)
print(f"\nThe largest immanant for this matrix is: {largest_imm}")

<<<0>>>