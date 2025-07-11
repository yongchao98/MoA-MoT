def solve_math_problem():
    """
    This function calculates the total mass based on the problem description.
    Since the variables q and q_v are symbolic, we will represent them as strings
    and print the derivation of the final formula.
    """
    # Symbolic representation of the variables
    q = "q"
    q_v = "q_v"

    # Step 1: State the volume of the measure mu.
    # The total mass of the measure mu on the space GL_n^1(K_hat) / GL_n(R),
    # denoted as vol(mu), is 1/(q-1). This result is derived from the theory of
    # S-arithmetic groups over function fields, leveraging the fact that the
    # Tamagawa number of SL_n is 1.
    volume_mu_str = f"1 / ({q} - 1)"

    # Step 2: State the factor given in the problem.
    factor_str = f"({q_v} * ({q} - 1)) / ({q_v} - 1)"

    # Step 3: Perform the final calculation and print the steps.
    print("The problem is to find the total mass of the measure resulting from the product:")
    print(f"({factor_str}) * mu")
    print("\nStep-by-step calculation:")
    print(f"1. The total mass of the base measure mu is: vol(mu) = {volume_mu_str}")
    print(f"2. The total mass is the product of the factor and vol(mu):")
    print(f"   Total Mass = ({factor_str}) * ({volume_mu_str})")
    print("\n3. We can write this as a single fraction:")
    print(f"   Total Mass = ({q_v} * ({q} - 1)) / (({q_v} - 1) * ({q} - 1))")
    print("\n4. The term '({q} - 1)' appears in both the numerator and the denominator, so it cancels out.")
    final_result_str = f"{q_v} / ({q_v} - 1)"
    print(f"5. The final simplified result is:")
    print(f"   Total Mass = {final_result_str}")

# Execute the function to print the solution.
solve_math_problem()