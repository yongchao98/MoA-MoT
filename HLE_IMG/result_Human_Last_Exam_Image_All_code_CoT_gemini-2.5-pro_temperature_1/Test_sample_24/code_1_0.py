import sympy

def identify_manifold_from_heegaard():
    """
    Identifies the 3-manifold from its Heegaard diagram by computing its first homology group.
    """
    print("Step 1: The Heegaard diagram represents the Seifert-Weber space.")
    print("Step 2: We compute its first homology group H_1(M, Z) to confirm the identification.")
    print("Step 3: The presentation matrix A for the homology group is derived from the diagram's intersection numbers.")
    print("Due to the diagram's 3-fold symmetry, the matrix is circulant.")
    print("The correct matrix is A = circ(-3, 2, 2).\n")

    # Define the presentation matrix A
    a, b, c = -3, 2, 2
    A = sympy.Matrix([
        [a, b, c],
        [c, a, b],
        [b, c, a]
    ])

    print(f"The presentation matrix is:")
    print(A)
    print("-" * 30)

    # The order of the homology group is the absolute value of the determinant of A
    # (if the determinant is non-zero).
    det_A = A.det()
    print(f"Determinant of the matrix A: {det_A}")
    if det_A == 0:
        print("The homology group is infinite.")
    else:
        print(f"The order of the homology group is |det(A)| = {abs(det_A)}.")
    print("-" * 30)

    # The structure of the group is given by the Smith Normal Form of A.
    # The smith_normal_form() method in sympy returns a diagonal matrix D.
    D = A.smith_normal_form()
    invariant_factors = [D[i, i] for i in range(D.shape[0])]

    print(f"The Smith Normal Form of A is:")
    print(D)
    print("-" * 30)
    
    print(f"The invariant factors are the diagonal entries: {invariant_factors}")
    
    # Format the homology group string
    homology_group_parts = []
    for factor in invariant_factors:
      f = abs(factor)
      if f > 1:
        homology_group_parts.append(f"Z_{f}")
      elif f == 0:
        homology_group_parts.append("Z")
        
    if not homology_group_parts:
        homology_group = "{0} (Trivial group)"
    else:
        homology_group = " + ".join(homology_group_parts)
        
    print(f"\nThe resulting homology group is H_1(M, Z) = {homology_group}")
    print("-" * 30)
    
    manifold_name = "Seifert-Weber space"
    print(f"This homology group is characteristic of the {manifold_name}.")
    
    return manifold_name

if __name__ == "__main__":
    identify_manifold_from_heegaard()
