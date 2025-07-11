def solve_controller():
    """
    This function determines and prints the set of all proper stabilizing controllers H_2(s)
    for the system H_1(s) = s / (s^2 - 1).
    The solution is derived using Youla-Kucera parametrization.
    """
    
    # The numerator of H_2(s) is parametrized by K(s)
    numerator = "4*s**2 + 8*s + 4 + (s**2 - 1)*K(s)"
    
    # The denominator of H_2(s) is parametrized by K(s)
    denominator = "s**2 - 1 - s*K(s)"
    
    # Print the resulting transfer function for H_2(s)
    # The prompt requires us to output each number in the final equation.
    print("The set of all proper stabilizing controllers is given by H_2(s):")
    print(f"H_2(s) = ({numerator}) / ({denominator})")
    print("\nwhere K(s) is any stable and proper rational function.")

solve_controller()