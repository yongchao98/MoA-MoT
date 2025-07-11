def solve_controller():
    """
    This function prints the transfer function for the set of all proper
    stabilizing controllers H_2(s) for the plant H_1(s) = s / (s^2 - 1).
    """

    # Numerator and Denominator of the plant H_1(s)
    N_s = "s"
    D_s = "s^2-1"

    # A suitable solution pair for the Bezout Identity N*X + D*Y = 1
    # X(s) = -s^2 + s + 1
    # Y(s) = s - 1
    X_s_coeff_2 = -1
    X_s_coeff_1 = 1
    X_s_coeff_0 = 1

    Y_s_coeff_1 = 1
    Y_s_coeff_0 = -1

    # The general form of the controller H_2(s) is (X + D*K) / (Y - N*K)
    # where K(s) is a stable and strictly proper rational function parameter.

    # We will print the equation for H_2(s) with all its coefficients.
    
    # Numerator of H_2(s)
    num_str = f"({X_s_coeff_2})s^2 + ({X_s_coeff_1})s + ({X_s_coeff_0}) + ({D_s})K(s)"
    
    # Denominator of H_2(s)
    den_str = f"({Y_s_coeff_1})s + ({Y_s_coeff_0}) - ({N_s})K(s)"
    
    print("The set of all proper stabilizing controllers is given by H_2(s):")
    print("")
    # Print the fraction line by line for clarity
    print(f"         {num_str}")
    print("H_2(s) = " + "-" * (len(num_str) + 2))
    print(f"         {den_str}")
    print("")
    print("For H_2(s) to be a proper controller, K(s) must be any stable and strictly proper rational function.")

solve_controller()