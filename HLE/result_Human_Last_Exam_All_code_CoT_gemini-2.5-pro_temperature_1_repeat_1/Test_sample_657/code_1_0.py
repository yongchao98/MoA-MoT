def display_controller_formula():
    """
    This function prints the set of all proper stabilizing controllers H_2(s)
    for the plant H_1(s) = s / (s^2 - 1).
    The controller is parametrized by a function K(s).
    """
    
    # Define the numerator and denominator of the controller H_2(s) as strings.
    # The controller is of the form: (B(s)K(s) + A(s)) / (D(s)K(s) + C(s))
    
    # Numerator: (s^2 - 1)K(s) + 4s^2 + 8s + 4
    numerator_str = "(s**2 - 1)*K(s) + 4*s**2 + 8*s + 4"
    
    # Denominator: -s*K(s) + s^2 - 1
    denominator_str = "-s*K(s) + s**2 - 1"
    
    print("The set of all proper stabilizing controllers H_2(s) is given by:")
    print("-" * 60)
    print("H_2(s) = ( {} ) / ( {} )".format(numerator_str, denominator_str))
    print("-" * 60)
    print("\nWhere K(s) is any stable and proper rational function.")
    print("This means:")
    print("1. All poles of K(s) must be in the left-half of the complex plane.")
    print("2. The degree of the numerator of K(s) must be less than or equal to the degree of its denominator.")

if __name__ == '__main__':
    display_controller_formula()