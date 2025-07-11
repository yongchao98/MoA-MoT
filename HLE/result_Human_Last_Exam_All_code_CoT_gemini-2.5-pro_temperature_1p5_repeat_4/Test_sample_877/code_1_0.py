import sympy

def solve():
    """
    This function derives and prints the function h(x).
    """
    x = sympy.Symbol('x')

    # The function h(x) is derived from the separatrix equation
    h_x = 4*x**2 + 2*x*sympy.ln(2*x) - 6*x + 2

    print("The function h(x) is given by the expression:")
    print(h_x)
    
    # The problem asks to output each number in the final equation.
    # We can represent the function as h(x) = c1*x**2 + c2*x*ln(2*x) + c3*x + c4
    c1 = 4
    c2 = 2
    c3 = -6
    c4 = 2
    
    print("\nIn a structured format, the equation is h(x) = a*x^2 + b*x*ln(2x) + c*x + d")
    print(f"The coefficients are:")
    print(f"a = {c1}")
    print(f"b = {c2}")
    print(f"c = {c3}")
    print(f"d = {c4}")

solve()