def solve_controller():
    """
    This function prints the set of all proper stabilizing controllers H_2(s)
    for the plant H_1(s) = s/(s^2-1).
    The result is parametrized by a function K(s), which must be any stable,
    proper rational function for the controller H_2(s) to be proper.
    """

    # The numerator of H_2(s) is parametrized by K(s)
    numerator = "4*s**2 + 8*s + 4 + (s**2 - 1)*K(s)"
    
    # The denominator of H_2(s) is parametrized by K(s)
    denominator = "s**2 - 1 - s*K(s)"
    
    # Printing each part of the final equation
    print("The set of all proper stabilizing controllers is given by H_2(s) = N_c(s) / D_c(s), where K(s) is any stable and proper rational function.")
    print("\nNumerator N_c(s):")
    print(numerator)
    print("\nDenominator D_c(s):")
    print(denominator)
    
    # Constructing and printing the final transfer function H_2(s)
    print("\nFull transfer function H_2(s):")
    
    # To make the output clear, we find the length of the longest string (num or den)
    # and use it to format the fraction bar nicely.
    max_len = max(len(numerator), len(denominator))
    fraction_bar = "-" * max_len
    
    print(f"H_2(s) = {numerator}")
    print(f"         {fraction_bar}")
    print(f"         {denominator}")

solve_controller()