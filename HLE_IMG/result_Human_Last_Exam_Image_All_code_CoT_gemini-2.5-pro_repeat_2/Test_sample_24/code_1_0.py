import sympy

def identify_manifold_from_heegaard():
    """
    Identifies a three-manifold from its Heegaard diagram by calculating
    its first homology group.
    """
    # Step 1: Define the presentation matrix from the relations derived from the diagram.
    # The generators are g1, g2, g3 corresponding to α1, α2, α3.
    # The relations are:
    # g1 + g2 = 0  (from β3)
    # g2 + g3 = 0  (from β1)
    # g1 + g3 = 0  (from β2)
    # We reorder them for clarity.
    relations = [
        (1, 1, 0),  # g1 + g2
        (0, 1, 1),  # g2 + g3
        (1, 0, 1)   # g1 + g3
    ]
    A = sympy.Matrix(relations)

    print("The presentation matrix for the first homology group H_1(M) is:")
    print(A)
    print("-" * 30)

    # Step 2: Compute the Smith Normal Form of the matrix.
    # The diagonal entries of the Smith Normal Form (d1, d2, ...) determine the
    # structure of the homology group as a direct sum of cyclic groups:
    # H_1(M, Z) = Z/d1Z + Z/d2Z + ...
    D = A.smith_normal_form()

    print("The Smith Normal Form of the matrix is:")
    print(D)
    print("-" * 30)

    # Step 3: Interpret the result to identify the homology group.
    diag_elements = [D[i, i] for i in range(D.shape[0])]

    print("The structure of the first homology group is determined by the diagonal elements.")
    # The final equation with all numbers:
    print(f"H_1(M, Z) = Z/{diag_elements[0]}Z + Z/{diag_elements[1]}Z + Z/{diag_elements[2]}Z")

    # Simplify the expression for the group
    homology_group_parts = []
    for d in diag_elements:
        d_int = abs(int(d))
        if d_int > 1:
            homology_group_parts.append(f"Z_{d_int}")
    
    if not homology_group_parts:
        homology_group = "0" # The trivial group
    else:
        homology_group = " + ".join(homology_group_parts)

    print(f"Thus, H_1(M, Z) is isomorphic to {homology_group}.")
    print("-" * 30)

    # Step 4: Identify the manifold.
    print("Conclusion:")
    print("The manifold is a closed, orientable 3-manifold with its first homology group being Z_2.")
    print("The simplest and most common manifold with this property is the Lens Space L(2,1),")
    print("which is also known as the real projective 3-space (RP^3).")

if __name__ == '__main__':
    identify_manifold_from_heegaard()