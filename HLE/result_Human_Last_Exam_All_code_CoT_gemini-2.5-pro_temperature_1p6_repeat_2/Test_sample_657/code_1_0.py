import sympy

def find_stabilizing_controllers():
    """
    Determines the set of all stabilizing controllers H_2(s) for the system H_1(s) = s/(s^2-1)
    using Youla-Kucera parametrization.
    """
    # Define the symbolic variable 's' and the parameter function 'K(s)'
    s = sympy.symbols('s')
    K = sympy.Function('K')(s)

    # Define the plant H_1(s) = N(s) / D(s)
    N = s
    D = s**2 - 1

    # Find a particular solution (X(s), Y(s)) to the Bezout identity:
    # N(s)X(s) + D(s)Y(s) = 1
    # s*X(s) + (s^2 - 1)*Y(s) = 1
    # A simple solution is X(s) = s and Y(s) = -1.
    # Check: s*(s) + (s^2 - 1)*(-1) = s^2 - s^2 + 1 = 1.
    X = s
    Y = -1

    # Youla-Kucera parametrization for the set of all stabilizing controllers H_2(s) is:
    # H_2(s) = (X(s) + D(s)K(s)) / (Y(s) - N(s)K(s))
    # where K(s) is any stable and proper rational function.
    
    # Construct the numerator and denominator of H_2(s)
    H2_numerator = X + D * K
    H2_denominator = Y - N * K
    
    # To make the denominator look a bit cleaner, we can factor out the minus sign
    H2_denominator_simplified = -( -Y + N * K)

    # Print the result in a readable format
    print("The system is given by H_1(s) = N(s)/D(s):")
    print(f"N(s) = {N}")
    print(f"D(s) = {D}")
    print("\nA solution to the Bezout identity N(s)X(s) + D(s)Y(s) = 1 is:")
    print(f"X(s) = {X}")
    print(f"Y(s) = {Y}")
    print("\nThe set of all stabilizing controllers H_2(s) is parametrized by a stable, proper function K(s):")
    print("H_2(s) = [X(s) + D(s)K(s)] / [Y(s) - N(s)K(s)]")
    
    # Display the final parametrized transfer function for H_2(s)
    # Using sympy.pretty_print for a nicer output
    print("\nSubstituting the expressions for N, D, X, Y gives:")
    
    final_num_expr = sympy.collect(H2_numerator, K)
    final_den_expr = sympy.collect(H2_denominator_simplified, K)
    
    # We use 'print' to show the step-by-step thinking for the final equation format.
    print(f"Numerator = {final_num_expr}")
    print(f"Denominator = {final_den_expr}")

    # Final answer display format, still outputting each part of the final equation.
    print("\nFinal Answer:")
    print("H_2(s) = ({}) / ({})".format(final_num_expr, H2_denominator))
    

find_stabilizing_controllers()