import sympy

def solve_controller():
    """
    This function uses sympy to define and print the set of all proper stabilizing
    controllers H_2(s) for the given system H_1(s).
    """
    # Define the symbolic variable 's' for the transfer function
    s = sympy.Symbol('s')

    # Define K as a symbolic function of s. In control theory, K(s) must be a
    # stable and proper rational function for H_2(s) to be a stabilizing controller.
    K = sympy.Function('K')(s)

    # Construct the numerator of the controller H_2(s) based on the Youla-Kucera parametrization
    # Numerator = X_0*(s+1)^2 + D*(s+1)^2*K = 4*(s+1)**2 + (s-1)*(s+1)*K
    num_H2 = 4*s**2 + 8*s + 4 + (s**2 - 1)*K

    # Construct the denominator of the controller H_2(s)
    # Denominator = Y_0*(s+1)^2 - N*(s+1)^2*K = (s-1)*(s+1) - s*K
    den_H2 = s**2 - 1 - s*K

    # Collect terms with K(s) for a structured view
    num_collected = sympy.collect(num_H2, K)
    den_collected = sympy.collect(den_H2, K)

    print("The set of all proper stabilizing controllers H_2(s) is given by the formula:")
    
    # Create a string representation for manual printing to ensure clarity
    num_str = f"({num_collected})"
    den_str = f"({den_collected})"

    # Print the final result in a clear, readable format
    print(f"\n          {num_str}")
    print(f"H_2(s) = {'-' * len(num_str)}")
    print(f"          {den_str}")
    
    print("\nWhere K(s) is any stable and proper rational function.")
    print("\nTo show all numbers in the equation, the expanded form is:")
    print(f"Numerator:   4*s**2 + 8*s + 4 + (s**2 - 1)*K(s)")
    print(f"Denominator: s**2 - 1 - s*K(s)")

solve_controller()