import sympy

def find_stabilizing_controllers():
    """
    This function derives and prints the parametrization of all stabilizing controllers
    for the plant H_1(s) = s / (s^2 - 1).
    """
    # Define the symbolic variable s and the parameter function K(s)
    s = sympy.Symbol('s')
    
    # The polynomials from the Bezout identity: N*X + D*Y = 1
    # For H_1(s) = s/(s^2-1), we have N=s, D=s^2-1, X=s, Y=-1.
    X = "s"
    Y = "-1"
    N = "s"
    D = "s^2 - 1"
    
    # Construct the numerator and denominator of the controller H_2(s)
    # H_2(s) = (X + D*K) / (Y - N*K)
    numerator_str = f"{X} + ({D})*K(s)"
    denominator_str = f"{Y} - {N}*K(s)"
    
    # Determine the length for the separator line for pretty printing
    separator_len = max(len(numerator_str), len(denominator_str))
    separator = "-" * separator_len
    
    print("The set of all stabilizing controllers H_2(s) is given by:")
    print()
    print(f"H_2(s) =  {numerator_str}")
    print(f"          {separator}")
    print(f"          {denominator_str}")
    print()
    print("where K(s) is any stable and proper rational function chosen")
    print("such that H_2(s) is also a proper transfer function.")

find_stabilizing_controllers()