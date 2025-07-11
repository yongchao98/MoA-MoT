import sympy

def solve_manifold_identity():
    """
    This function analyzes the Heegaard diagram to identify the 3-manifold.
    It constructs the intersection matrix and computes its Smith Normal Form
    to determine the first homology group of the manifold.
    """
    # Step 1: Define the intersection matrix based on the diagram analysis.
    # A[i, j] = algebraic intersection number of beta_i and alpha_j.
    # beta_1 ~ alpha_2 - alpha_3
    # beta_2 ~ alpha_3 - alpha_1
    # beta_3 ~ alpha_1 - alpha_2
    A = sympy.Matrix([
        [0, 1, -1],
        [-1, 0, 1],
        [1, -1, 0]
    ])

    print("Step 1: The Heegaard diagram analysis leads to the following intersection matrix A:")
    print("A =")
    sympy.pprint(A)
    print("\nEach row represents a beta-curve's homology class in terms of the alpha-curves' classes.")
    print(f"For example, row 1 [0, 1, -1] corresponds to the relation [β₁] = 0*α₁ + 1*α₂ - 1*α₃.\n")


    # Step 2: Compute the Smith Normal Form of the matrix.
    # The Smith Normal Form diagonalizes the matrix over integers.
    # The diagonal entries determine the structure of the homology group.
    print("Step 2: Compute the Smith Normal Form (SNF) of matrix A.")
    snf_matrix = A.smith_normal_form()
    print("The Smith Normal Form of A is:")
    sympy.pprint(snf_matrix)
    print("")

    # Step 3: Interpret the result to find the first homology group H_1(M).
    # H_1(M) = Z/d1 Z + Z/d2 Z + ...
    # where di are the diagonal elements of the SNF.
    # Z/0Z = Z and Z/1Z = {0} (the trivial group).
    diag_elems = [snf_matrix[i, i] for i in range(snf_matrix.rows)]
    homology_group = []
    for d in diag_elems:
        if d == 0:
            homology_group.append("Z")
        elif d != 1:
            homology_group.append(f"Z/{d}Z")

    if not homology_group:
        homology_group.append("{0}")

    homology_str = " \u2295 ".join(homology_group) # \u2295 is the direct sum symbol

    print("Step 3: Determine the first homology group H₁(M, Z) from the SNF.")
    print(f"The diagonal elements are {diag_elems[0]}, {diag_elems[1]}, and {diag_elems[2]}.")
    print("The homology group H₁(M) is the direct sum of cyclic groups Z/d_i*Z.")
    print(f"H₁(M, Z) = Z/{diag_elems[0]}Z \u2295 Z/{diag_elems[1]}Z \u2295 Z/{diag_elems[2]}Z")
    print(f"Since Z/1Z is the trivial group {0} and Z/0Z is Z, we get:")
    print(f"H₁(M, Z) = {homology_str}\n")

    # Step 4: Identify the manifold.
    print("Step 4: Identify the manifold.")
    print("The manifold M is a compact, orientable 3-manifold with H₁(M, Z) = Z.")
    print("The simplest such manifold is the product of a circle and a 2-sphere, S¹ \u00d7 S².")
    print("\nConclusion: The Heegaard diagram represents the 3-manifold S¹ \u00d7 S².")

solve_manifold_identity()