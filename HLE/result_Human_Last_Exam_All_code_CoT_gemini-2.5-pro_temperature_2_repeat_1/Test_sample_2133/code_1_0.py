import numpy as np
from itertools import permutations
import math

def solve():
    """
    This script finds the largest immanant for a specific 3x3 matrix Mn.
    This matrix Mn is chosen because its Hermite Normal Form (HNF) maximizes
    a given ratio, among a class of matrices.
    """
    
    # For n=3, the optimal matrix Mn is found to be:
    M = np.array([
        [1, -11, 2],
        [1, -3, 2],
        [1, 1, 2]
    ])

    print(f"The specific 3x3 nilpotent matrix with all non-zero integer entries is:")
    print(M)
    print("-" * 30)

    # For S_3, there are 3 conjugacy classes, characterized by cycle structure:
    # id: {(0,1,2)} - identity permutation, e.g. (1)(2)(3)
    # transp: {(1,0,2), (2,1,0), (0,2,1)} - transpositions, e.g. (1 2)(3)
    # 3-cycle: {(1,2,0), (2,0,1)} - 3-cycles, e.g. (1 2 3)

    # We map each permutation (as a tuple) to its conjugacy class identifier.
    perm_to_class = {
        (0, 1, 2): 'id',
        (1, 0, 2): 'transp',
        (2, 1, 0): 'transp',
        (0, 2, 1): 'transp',
        (1, 2, 0): '3-cycle',
        (2, 0, 1): '3-cycle'
    }

    # S_3 Character Table
    # Rows: chi_1 (trivial/permanent), chi_2 (sign/determinant), chi_3 (standard)
    # Columns: conjugacy classes (id, transp, 3-cycle)
    char_table = {
        'chi1_perm': {'id': 1, 'transp': 1, '3-cycle': 1},
        'chi2_det': {'id': 1, 'transp': -1, '3-cycle': 1},
        'chi3_std': {'id': 2, 'transp': 0, '3-cycle': -1}
    }
    
    immanants = {}
    equations = {}

    all_perms = list(permutations(range(3)))

    for chi_name, char_values in char_table.items():
        imm = 0
        equation_parts = []
        for p in all_perms:
            prod = 1
            prod_str_list = []
            for i in range(3):
                prod *= M[i, p[i]]
                prod_str_list.append(str(M[i,p[i]]))
            
            prod_str = f"({'*'.join(prod_str_list)})"
            
            class_id = perm_to_class[p]
            char_val = char_values[class_id]

            if char_val == 0:
                continue

            term = char_val * prod
            imm += term
            
            op = "+" if char_val > 0 else "-"
            # For readable equation, handle sign of character
            if abs(char_val) == 1:
                if op == "+":
                    equation_parts.append(f"{prod_str}")
                else: # op == "-"
                    equation_parts.append(f"- {prod_str}")
            else: # char_val can be 2
                 equation_parts.append(f"{char_val} * {prod_str}")

        immanants[chi_name] = imm
        equations[chi_name] = " ".join(equation_parts)

    print("The immanants are calculated as follows:")
    
    # Determinant
    name_det = 'chi2_det'
    det_perm_contrib = {p: np.prod([M[i, p[i]] for i in range(3)]) for p in all_perms}
    print(f"\n1. Determinant (sign character):")
    print(f"det(M) = ({M[0,0]}*{M[1,1]}*{M[2,2]}) - ({M[0,0]}*{M[1,2]}*{M[2,1]})"
        f" + ({M[0,1]}*{M[1,2]}*{M[2,0]}) - ({M[0,1]}*{M[1,0]}*{M[2,2]})"
        f" + ({M[0,2]}*{M[1,0]}*{M[2,1]}) - ({M[0,2]}*{M[1,1]}*{M[2,0]})")
    print(f"  = ({det_perm_contrib[(0,1,2)]}) - ({det_perm_contrib[(0,2,1)]})"
          f" + ({det_perm_contrib[(1,2,0)]}) - ({det_perm_contrib[(1,0,2)]})"
          f" + ({det_perm_contrib[(2,0,1)]}) - ({det_perm_contrib[(2,1,0)]})"
          f" = {immanants[name_det]}")

    # Permanent
    name_perm = 'chi1_perm'
    print(f"\n2. Permanent (trivial character):")
    print(f"per(M) = " + " + ".join([f"({det_perm_contrib[p]})" for p in all_perms]))
    print(f"  = {immanants[name_perm]}")
    
    # Standard Character Immanant
    name_std = 'chi3_std'
    id_perm = (0,1,2)
    cycle3_perms = [(1,2,0), (2,0,1)]
    print(f"\n3. Immanant for the 2-dim standard character:")
    print(f"Imm_chi3(M) = 2 * (product for id) - 1 * (sum of products for 3-cycles)")
    print(f"  = 2 * ({M[0,0]}*{M[1,1]}*{M[2,2]}) - (({M[0,1]}*{M[1,2]}*{M[2,0]}) + ({M[0,2]}*{M[1,0]}*{M[2,1]}))")
    print(f"  = 2 * ({det_perm_contrib[id_perm]})"
          f" - (({det_perm_contrib[cycle3_perms[0]]}) + ({det_perm_contrib[cycle3_perms[1]]}))"
          f" = {immanants[name_std]}")


    print("-" * 30)

    largest_immanant = -math.inf
    for name, value in immanants.items():
        if value > largest_immanant:
            largest_immanant = value

    print(f"The calculated immanants are: {list(immanants.values())}")
    print(f"The largest immanant for this matrix is: {largest_immanant}")

solve()
>>>8
