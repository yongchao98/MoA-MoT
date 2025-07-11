def solve_controller():
    """
    This function prints the formula for the set of all stabilizing controllers H_2(s)
    for the plant H_1(s) = s / (s^2 - 1).
    The controller is parametrized by a function K(s).
    """
    
    # Numerator of the controller H_2(s)
    # N_H2(s) = X(s) + D(s)K(s)
    # X(s) = s
    # D(s) = s^2 - 1
    num_str = "s + (s**2 - 1)*K(s)"
    
    # Denominator of the controller H_2(s)
    # D_H2(s) = Y(s) - N(s)K(s)
    # Y(s) = -1
    # N(s) = s
    den_str = "-1 - s*K(s)"
    
    # Print the final transfer function H_2(s)
    print("H_2(s) = ({}) / ({})".format(num_str, den_str))

solve_controller()