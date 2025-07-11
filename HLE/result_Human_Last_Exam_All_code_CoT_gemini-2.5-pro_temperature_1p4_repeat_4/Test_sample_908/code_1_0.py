import sympy

def calculate_charges():
    """
    This function calculates and prints the symbolic expressions for the total volume
    and surface charges based on the formulas provided in the selected answer choice.
    """
    # Define the symbolic variables
    V = sympy.Symbol('V')
    epsilon = sympy.Symbol('varepsilon')
    pi = sympy.pi
    L = sympy.Symbol('L')
    a = sympy.Symbol('a')
    b = sympy.Symbol('b')

    # Formulas from Answer Choice B
    # Total volume charge q_v
    q_v = -4 * V * epsilon * pi * L / (1 - a**2/b**2)
    
    # Total surface charge on the inner electrode q_s(a)
    q_s_a = 2 * pi * L * V * epsilon / (1 - a**2/b**2)

    # Total surface charge on the outer electrode q_s(b)
    q_s_b = -4 * pi * L * V * epsilon * a**2 / (b**2 * (1 - a**2/b**2))

    # Print the results in a formatted way
    print("Based on the provided answer choice B, the charges are:")
    
    print("\nTotal volume charge (q_v):")
    # To display it in a cleaner way, we can simplify the fraction
    q_v_simplified = -4 * V * epsilon * pi * L * b**2 / (b**2 - a**2)
    print(f"q_v = {sympy.pretty(q_v)}")
    
    print("\nTotal surface charge on inner electrode (q_s(r = a)):")
    q_s_a_simplified = 2 * pi * L * V * epsilon * b**2 / (b**2 - a**2)
    print(f"q_s(a) = {sympy.pretty(q_s_a)}")

    print("\nTotal surface charge on outer electrode (q_s(r = b)):")
    q_s_b_simplified = -4 * pi * L * V * epsilon * a**2 / (b**2 - a**2)
    print(f"q_s(b) = {sympy.pretty(q_s_b)}")

if __name__ == '__main__':
    calculate_charges()