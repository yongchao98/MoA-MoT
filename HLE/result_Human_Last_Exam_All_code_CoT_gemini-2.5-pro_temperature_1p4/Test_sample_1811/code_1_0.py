import sympy

def solve_manifold_zeros():
    """
    Determines and explains the minimum number of zeros for a vector field on a
    compact manifold M with a non-empty boundary ∂M.
    """

    # Define symbolic variables for the Euler characteristics of the manifold M
    # and its boundary ∂M. These are integers.
    chi_M = sympy.Symbol('χ(M)')
    chi_dM = sympy.Symbol('χ(∂M)')

    # Step 1: State the condition for the existence of a non-vanishing vector field.
    # A non-vanishing vector field (one with 0 zeros) can exist on M if and only if
    # the Euler characteristic of its double, D(M), is zero.
    print("The minimum number of zeros depends on a topological condition.")
    print("A vector field with 0 zeros exists if and only if the Euler characteristic of the manifold's double, χ(D(M)), is 0.")

    # Step 2: Express the Euler characteristic of the double in terms of χ(M) and χ(∂M).
    # The formula is χ(D(M)) = 2*χ(M) - χ(∂M).
    print("\nThe condition is determined by the following equation:")
    print("Condition C = 2 * χ(M) - χ(∂M)")
    
    # Printing the numbers in the equation as requested
    number_2 = 2
    
    print(f"The numbers in this equation are: {number_2}")
    print("\nIf C = 0, the minimum number of zeros is 0.")

    # Step 3: Address the case where a non-vanishing field does not exist.
    # If C ≠ 0, any vector field must have at least one zero.
    # A theorem by Hopf states that a vector field with a single zero can be constructed.
    print("If C ≠ 0, any vector field must have at least one zero.")
    print("In this case, it is always possible to construct a vector field with exactly one zero.")
    
    number_1 = 1
    number_0 = 0

    print(f"Therefore, the minimum number of zeros is {number_1}.")

    # Step 4: Summarize the final answer.
    print("\n---")
    print("Final Answer Summary:")
    print("The least number of zeros a vector field can have on M is:")
    print(f"- {number_0}, if 2 * χ(M) - χ(∂M) = {number_0}")
    print(f"- {number_1}, if 2 * χ(M) - χ(∂M) ≠ {number_0}")

solve_manifold_zeros()