import sympy
from sympy.combinatorics.permutations import Permutation
from sympy.combinatorics.perm_groups import SymmetricGroup

def solve_matrix_problem():
    """
    Solves the problem for n=3 by analyzing a specific Mercer matrix.
    
    The chosen matrix M3 is constructed by a similarity transformation of the
    3x3 nilpotent Jordan block, ensuring it is 3-nilpotent and has
    non-zero integer entries. The function then calculates all immanants
    of this matrix and finds the largest one by magnitude.
    """
    n = 3

    # Define the specific Mercer matrix for n=3
    # This matrix M is 3-nilpotent and has all non-zero integer entries.
    # M = S * J * S^-1 where J is the Jordan block and S is a unimodular matrix.
    # S = Matrix([[1, 1, 1], [1, 2, 1], [1, 1, 2]])
    # J = Matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]])
    # S_inv = S.inv()
    # M_calc = S * J * S_inv
    # This yields M_calc = Matrix([[-2, 1, 1], [-3, 1, 2], [-2, 1, 1]])
    M = sympy.Matrix([[-2, 1, 1],
                      [-3, 1, 2],
                      [-2, 1, 1]])

    # Get the symmetric group S_n
    G = SymmetricGroup(n)
    
    # Get the character table for S_n
    # The columns correspond to the conjugacy classes and rows to irreducible characters.
    char_table = G.character_table()
    
    # Store all permutation objects, grouped by conjugacy class
    conjugacy_classes_perms = [list(cc) for cc in G.conjugacy_classes()]

    # We will compute each term of the immanant formula: prod(M[i, sigma(i)])
    # and store it in a dictionary mapping permutation to the product term.
    perm_product_map = {}
    for perm_obj in G.elements:
        # Permutation objects in Sympy might be 0-indexed internally
        # but apply to indices 0 to n-1. Let's be careful.
        # The permutation needs to map row i to col sigma(i).
        # In SymPy, perm_obj(i) gives the element at pos i moved to.
        # So we need to compute the permutation as a list of images.
        perm_list = perm_obj.array_form
        prod = 1
        for i in range(n):
            prod *= M[i, perm_list[i]]
        perm_product_map[perm_obj] = prod

    print(f"The chosen Mercer matrix M for n={n} is:")
    print(M)
    print("\n" + "="*30 + "\n")
    print("Calculating immanants for each irreducible character of S3:")
    
    immanants = []
    
    # Iterate through each character (each row of the character table)
    for i in range(len(char_table.table)):
        char = char_table.table[i] # This is a Character object
        imm = 0
        # Sum over all conjugacy classes
        for j, cc_perms in enumerate(conjugacy_classes_perms):
            char_value = char.values[j]
            # Sum over all permutations in the current conjugacy class
            for perm_obj in cc_perms:
                 imm += char_value * perm_product_map[perm_obj]
        
        immanants.append(imm)
        char_name = f"chi_{i+1}"
        if i == 0:
            char_name += " (permanent)"
        elif char.values == sympy.polys.rootoftools.RootOf(sympy.Symbol('x')**2 - 1, 1) ** G.sgn(list(G.conjugacy_classes()[i])[0]):
             # This condition is complicated. A simpler check for sgn char
             if all(val in [-1, 1] for val in char.values):
                 char_name += " (determinant)"
        print(f"Immanant for character {char_name}: {imm}")

    # Find the largest immanant (by magnitude)
    if not immanants:
        largest_imm = "Not found"
    else:
        largest_imm = max(immanants, key=abs)
    
    print("\n" + "="*30 + "\n")
    print(f"The largest immanant (by magnitude) is: {largest_imm}")


solve_matrix_problem()