import sympy

def display_H_bound_formula():
    """
    This function derives and displays the formula for the upper bound H
    based on the analysis of the problem.
    """
    # Define the symbolic variables based on the problem description
    # H(a, b, c, d, r, t)
    # a corresponds to k
    # b corresponds to the L1 norm of rho at t=0
    # c corresponds to pi
    # d corresponds to nu
    # r corresponds to the positive lower bound of rho(tau, x)
    # t corresponds to time
    k, b, c, d, r, t = sympy.symbols('k b c d r t')

    # The problem states k < 0, so -k is positive.
    # The derived formula for the upper bound H is:
    # H = (-k * b * t) / (c * d**2 * r)
    H_formula = (-k * b * t) / (c * d**2 * r)

    # Print the result in a formatted way
    print("Based on the problem analysis and the assumption that 'r' is a positive lower bound for rho(tau, x),")
    print("the explicit upper bound H is given by the formula:")
    
    # We use unicode for nicer printing of mathematical symbols
    # pi -> π, nu -> ν, rho -> ρ
    symbol_map = {
        'k': 'k',
        'b': '||\u03C1(0,\u00B7)||_L1',
        'c': '\u03C0',
        'd': '\u03BD',
        'r': 'r',
        't': 't'
    }
    
    numerator = f"(-{symbol_map['k']}) * {symbol_map['b']} * {symbol_map['t']}"
    denominator = f"{symbol_map['c']} * {symbol_map['d']}\u00b2 * {symbol_map['r']}"

    # To create a fraction-like display
    width = max(len(numerator), len(denominator))
    print("\n" + numerator.center(width))
    print("H = " + "\u2014" * width)
    print(denominator.center(width) + "\n")
    
    print(f"Here, (-{symbol_map['k']}) is positive since {symbol_map['k']} < 0.")

display_H_bound_formula()