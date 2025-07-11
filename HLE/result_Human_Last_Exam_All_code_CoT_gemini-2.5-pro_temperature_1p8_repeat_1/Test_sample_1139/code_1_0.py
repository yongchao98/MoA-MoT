def solve_class_d_variables():
    """
    Calculates the number of non-Grassmann variables for the SUSY sigma-model
    of symmetry class D with a given number of replicas.
    """
    # Number of replicas as specified in the problem
    n = 2

    # The number of non-Grassmann (bosonic) variables for symmetry class D
    # is given by the dimension of the symmetric space O(4n)/(O(2n) x O(2n)).
    # This dimension simplifies to the formula: 4 * n^2.

    print("Problem: Find the number of non-Grassmann variables for the supersymmetric sigma-model.")
    print("Details:")
    print(f"  - Symmetry Class: D")
    print(f"  - Number of Replicas (n): {n}")
    print("-" * 40)
    print("The formula for the number of variables is dim(O(4n)/(O(2n) x O(2n))) = 4 * n^2.")
    print("\nCalculation steps:")

    # Perform the calculation
    result = 4 * (n ** 2)

    # Print the final equation with the numbers substituted, as requested.
    print(f"Number of variables = 4 * ({n}**2)")
    print(f"                    = 4 * {n**2}")
    print(f"                    = {result}")

    print("-" * 40)
    print(f"The final answer is: {result}")

solve_class_d_variables()