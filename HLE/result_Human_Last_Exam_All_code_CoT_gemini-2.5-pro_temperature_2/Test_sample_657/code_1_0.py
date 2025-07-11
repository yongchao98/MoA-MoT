import sympy as sp

def main():
    """
    This script calculates and displays the set of all proper stabilizing controllers
    for the plant H1(s) = s / (s^2 - 1).
    """
    
    # Define the symbolic variable 's' for the transfer function
    s = sp.symbols('s')

    # Define K(s), the free parameter for the Youla-Kucera parametrization.
    # K(s) can be any stable and proper transfer function.
    K = sp.Function('K')(s)

    # The controller H_2(s) is given by the formula:
    # H_2(s) = (X(s) + D(s)*K(s)) / (Y(s) - N(s)*K(s))
    # After calculation, this simplifies to the form:
    # H_2(s) = (num_poly_B(s) + num_poly_A(s)*K(s)) / (den_poly_D(s) + den_poly_C(s)*K(s))
    # where the polynomials are:
    num_poly_A = s**2 - 1
    num_poly_B = 8*s**2 + 41*s + 32
    den_poly_C = -s
    den_poly_D = s**2 - 16
    
    # Construct the numerator and denominator of the controller H_2(s)
    h2_numerator = num_poly_B + num_poly_A * K
    h2_denominator = den_poly_D + den_poly_C * K

    print("The set of all proper stabilizing controllers H_2(s) is given by:")
    
    # Use SymPy's pretty print for a clear mathematical representation
    # sp.init_printing(use_unicode=True)
    
    # To ensure clear output in all terminals, we format it manually.
    h2_expression = sp.Eq(sp.Function('H_2')(s), h2_numerator / h2_denominator, evaluate=False)
    
    # The following print statements show the final equation and all its numeric coefficients.
    print("\n H_2(s) = ( A(s) * K(s) + B(s) ) / ( C(s) * K(s) + D(s) )\n")
    print("Where K(s) is any stable and proper transfer function, and:")
    print(f"A(s) = {num_poly_A}")
    print(f"B(s) = {num_poly_B}")
    print(f"C(s) = {den_poly_C}")
    print(f"D(s) = {den_poly_D}")

    print("\nResulting in the final expression:")
    
    # Create an appealing string representation for the final answer
    final_num = f"({sp.pretty(num_poly_A, use_unicode=False)})*K(s) + ({sp.pretty(num_poly_B, use_unicode=False)})"
    final_den = f"({sp.pretty(den_poly_D, use_unicode=False)}) + ({sp.pretty(den_poly_C, use_unicode=False)})*K(s)"
    
    print(f" H_2(s) =  {final_num}")
    print("          " + "-" * (len(final_num)))
    print(f"           {final_den}")

if __name__ == '__main__':
    main()