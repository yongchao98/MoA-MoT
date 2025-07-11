def solve_knot_problem():
    """
    This function solves the problem using the sagemath library.
    Note: This code must be run in a SageMath environment.
    """
    try:
        from sage.all import BraidGroup, Knot, var
    except ImportError:
        print("This script requires the SageMath library. Please run it in a SageMath environment.")
        return

    # Define the variable for the polynomial
    z = var('z')

    # Step 1 & 2: Define the braid from the problem and create the knot from its closure
    print("--- Analyzing the first knot from the braid closure ---")
    
    # Braid group on 5 strands
    B5 = BraidGroup(5)
    
    # The braid word for beta is s_4^-1 s_4^-1 s_3^-1 s_4 s_3^-1 s_2 s_1^-1 s_3^-1 s_2^-1 s_2^-1 s_2^-1 s_1^-1
    # In numeric representation (s_i is i, s_i^-1 is -i):
    beta_word = [-4, -4, -3, 4, -3, 2, -1, -3, -2, -2, -2, -1]
    beta = B5(beta_word)
    
    # Create the knot K_beta from the closure of beta
    K_beta = Knot(beta)
    print(f"The braid beta corresponds to the knot: {K_beta.name()}")

    # Step 3: Calculate the Alexander-Conway polynomial for K_beta
    nabla_beta = K_beta.alexander_conway_polynomial()
    print(f"The Alexander-Conway polynomial for this knot is: \u2207(z) = {nabla_beta}")

    # Step 6 (Part 1): Extract the z^2 coefficient
    coeff_beta = nabla_beta.coefficient(z**2)
    print(f"The coefficient of z^2 is: {coeff_beta}")
    
    print("\n--- Analyzing the second knot, 10_4 ---")
    
    # Step 4: Create the knot 10_4
    K_10_4 = Knot('10_4')
    print(f"The second knot is: {K_10_4.name()}")
    
    # Step 5: Calculate the Alexander-Conway polynomial for 10_4
    nabla_10_4 = K_10_4.alexander_conway_polynomial()
    print(f"The Alexander-Conway polynomial for 10_4 is: \u2207(z) = {nabla_10_4}")

    # Step 6 (Part 2): Extract the z^2 coefficient
    coeff_10_4 = nabla_10_4.coefficient(z**2)
    print(f"The coefficient of z^2 is: {coeff_10_4}")

    # Step 7: Calculate and print the final difference
    print("\n--- Final Calculation ---")
    difference = coeff_beta - coeff_10_4
    print(f"The difference between the z^2 coefficients is: ({coeff_beta}) - ({coeff_10_4}) = {difference}")

solve_knot_problem()