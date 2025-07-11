import sympy
from sympy.combinatorics.permutations import Permutation
from sympy.combinatorics.perm_groups import SymmetricGroup

def solve_matrix_problem(n):
    """
    This function solves the problem for a given integer n.

    1. It constructs a specific n x n n-nilpotent matrix M_n with non-zero integer entries.
       This matrix is a candidate for maximizing the specified ratio.
       The matrix is constructed as M_n = S_n * J_n(0) * S_n^(-1), where J_n(0) is the
       n x n Jordan block for eigenvalue 0, and S_n is a specially chosen
       unimodular integer matrix.

    2. It then calculates all the immanants of this matrix M_n.

    3. Finally, it returns the largest immanant.
    """
    if n < 2:
        raise ValueError("n must be an integer greater than or equal to 2.")

    # Step 1: Construct the matrix M_n
    
    # Construct the unimodular matrix S_n
    S_list = [[1] * n for _ in range(n)]
    for i in range(1, n):
        for j in range(1, n):
            S_list[i][j] = 1
    for i in range(1, n):
        S_list[i][i] = 2
    S = sympy.Matrix(S_list)

    # Construct the Jordan block J_n(0)
    J = sympy.diag(1, 1, k=1)
    
    # Calculate M_n = S * J * S^(-1)
    S_inv = S.inv()
    M = S * J * S_inv

    print(f"For n = {n}, the constructed matrix M_n is:")
    sympy.pprint(M)
    print("-" * 30)

    # Step 2: Calculate all immanants of M_n
    
    # Get the symmetric group S_n and its character table
    G = SymmetricGroup(n)
    try:
        # sympy > 1.10
        char_table = G.character_table()
    except AttributeError:
        # older sympy
        char_table = G.get_character_table()


    # Get all permutations and their conjugacy classes
    perms = list(G.elements)
    conjugacy_classes = G.conjugacy_classes()
    
    # Map each permutation to its character values for all irreducible characters
    perm_to_chars = {}
    for p in perms:
        # Find which conjugacy class the permutation belongs to
        for i, cc in enumerate(conjugacy_classes):
            if p in cc:
                # Store the column of character values for this class
                perm_to_chars[p] = [row[i] for row in char_table]
                break

    immanants = []
    num_chars = len(char_table)

    for char_index in range(num_chars):
        imm = 0
        for p in perms:
            # Get the character value for the current character and permutation
            char_val = perm_to_chars[p][char_index]
            
            # Calculate the product term
            prod = 1
            perm_list = p.array_form
            for i in range(n):
                prod *= M[i, perm_list[i]]
            
            imm += char_val * prod
        immanants.append(imm)

    print(f"The immanants of M_{n} are:")
    print(immanants)
    print("-" * 30)

    # Step 3: Find and print the largest immanant
    largest_immanant = max(immanants)
    print(f"The largest immanant of M_{n} is: {largest_immanant}")
    return largest_immanant

if __name__ == '__main__':
    # The user did not specify n, so we choose n=4 as a representative example.
    n_val = 4
    final_answer = solve_matrix_problem(n_val)
    # The final answer format requested by the user.
    # print(f'<<<{final_answer}>>>')