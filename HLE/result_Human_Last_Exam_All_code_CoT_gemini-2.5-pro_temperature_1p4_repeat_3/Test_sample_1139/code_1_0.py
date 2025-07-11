def solve_class_d_variables():
    """
    Calculates the number of non-Grassmann variables for the SUSY sigma-model
    of symmetry class D with a given number of replicas.
    """
    # The user specified two replicas.
    n = 2

    # Step 1: Calculate the dimensions of the relevant Lie groups for n=2.
    # The groups are Sp(2n, R), O(2n), and U(n).
    dim_Sp2n_R = n * (2 * n + 1)
    dim_O2n = n * (2 * n - 1)
    dim_Un = n * n

    # Step 2: The bosonic manifold is M = (Sp(2n, R) / U(n)) x (O(2n) / U(n)).
    # We calculate the dimension of the non-compact and compact parts.

    # Non-compact part: Sp(2n, R) / U(n)
    dim_non_compact = dim_Sp2n_R - dim_Un

    # Compact part: O(2n) / U(n)
    dim_compact = dim_O2n - dim_Un

    # Step 3: The total number of variables is the sum of the dimensions
    # of the two parts of the manifold.
    total_dim = dim_non_compact + dim_compact

    # Print the explanation and the final calculation.
    print("The number of non-Grassmann variables is the dimension of the bosonic target manifold.")
    print(f"For symmetry class D with n={n} replicas, this manifold is M = (Sp({2*n}, R) / U({n})) x (O({2*n}) / U({n})).\n")
    print("The total dimension is the sum of the dimensions of its non-compact and compact parts.\n")

    print("Dimension of the non-compact part:")
    print(f"dim(Sp({2*n}, R) / U({n})) = dim(Sp({2*n}, R)) - dim(U({n}))")
    print(f"                     = {dim_Sp2n_R} - {dim_Un} = {dim_non_compact}\n")

    print("Dimension of the compact part:")
    print(f"dim(O({2*n}) / U({n})) = dim(O({2*n})) - dim(U({n}))")
    print(f"                   = {dim_O2n} - {dim_Un} = {dim_compact}\n")

    print("Total number of non-Grassmann variables:")
    print(f"Total Dim = (Dimension of non-compact part) + (Dimension of compact part)")
    print(f"Total Dim = {dim_non_compact} + {dim_compact} = {total_dim}")


solve_class_d_variables()