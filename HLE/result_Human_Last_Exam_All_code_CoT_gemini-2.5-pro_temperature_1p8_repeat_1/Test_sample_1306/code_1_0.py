def solve_quantum_representation_percentage():
    """
    Calculates the percentage of irreducible representations among all indecomposable
    representations for the Taft algebra T_l, which is a finite-dimensional
    quantum group analogue for sl_2.

    The problem specifies q as a primitive third root of unity, so l=3.
    """

    # For a primitive l-th root of unity, we have l=3.
    l = 3

    # The number of isomorphism classes of irreducible (simple) representations
    # for the Taft algebra T_l is l.
    num_irreducible = l

    # The number of isomorphism classes of indecomposable representations
    # for the Taft algebra T_l is l^2. This algebra is of finite representation type.
    num_indecomposable = l**2

    # Calculate the percentage.
    percentage = (num_irreducible / num_indecomposable) * 100

    print("This problem is interpreted in the context of the Taft algebra T_l, for which the number of indecomposable representations is finite.")
    print(f"For q a primitive 3rd root of unity, we have l = {l}.")
    print(f"The number of irreducible representations is equal to l, which is {num_irreducible}.")
    print(f"The number of indecomposable representations is equal to l^2, which is {num_indecomposable}.")
    print("\nThe percentage of irreducible representations is calculated as:")
    print(f"Percentage = 100 * (Number of Irreducible / Number of Indecomposable)")
    print(f"Percentage = 100 * ({num_irreducible} / {num_indecomposable}) = {percentage:.2f}%")

solve_quantum_representation_percentage()
