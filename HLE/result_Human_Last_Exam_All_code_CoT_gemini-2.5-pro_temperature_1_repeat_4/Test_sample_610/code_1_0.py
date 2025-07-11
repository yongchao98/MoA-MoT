import sympy

def solve():
    """
    Calculates the exact value of l(n, b) based on the derived formula.
    The user will be prompted to enter values for n and b.
    """
    # Using sympy for symbolic manipulation to keep the expression exact if possible.
    n_str = input("Enter the value of n (integer, n >= 10): ")
    b_str = input("Enter the value of b (float, 0 < b < 1): ")

    try:
        n = sympy.Integer(n_str)
        # Try to represent b as a rational number for exactness
        b = sympy.Rational(b_str)
    except (ValueError, TypeError):
        print("Invalid input. Please enter valid numerical values.")
        return

    if not n.is_Integer or n < 10:
        print("n must be an integer greater than or equal to 10.")
        return
    if not (b > 0 and b < 1):
        print("b must be in the interval (0, 1).")
        return

    # Derived formula for l(n, b)
    # l(n,b) = (2*(n-1) - 2*b*(n-2)) / (1+b)
    
    numerator = 2 * (n - 1) - 2 * b * (n - 2)
    denominator = 1 + b
    
    l_value = numerator / denominator

    # Simplify the expression
    simplified_l = sympy.simplify(l_value)

    print(f"The exact value of l(n, b) for n={n} and b={b} is:")
    
    # To show the equation for the final answer
    num_val = 2 * (n - 1)
    num_b_val = -2 * (n - 2)
    den_val = 1
    den_b_val = 1

    # Print the equation form
    # Using Sympy's pretty print for a clear mathematical representation
    n_sym, b_sym = sympy.symbols('n b')
    
    # We print the components of the formula as requested
    print(f"l(n,b) = (A + B*b) / (C + D*b)")
    print(f"A = {2*(n-1)}")
    print(f"B = {-2*(n-2)}")
    print(f"C = 1")
    print(f"D = 1")
    
    print("\nSubstituting the values n and b:")
    
    final_numerator = 2*(n-1) - 2*b*(n-2)
    final_denominator = 1+b
    
    print(f"l({n}, {b}) = ({final_numerator}) / ({final_denominator})")
    print(f"l({n}, {b}) = {simplified_l}")

    # For the final answer format
    # Since the problem is symbolic, we express the formula itself as the answer.
    # The format <<<answer>>> is for a single value, so we will output the formula structure.
    # The problem asks for the *exact value*, which in this symbolic case is the formula.
    # As a final symbolic representation:
    # l(n,b) = (2n - 2 - (2n-4)b) / (1+b)

solve()