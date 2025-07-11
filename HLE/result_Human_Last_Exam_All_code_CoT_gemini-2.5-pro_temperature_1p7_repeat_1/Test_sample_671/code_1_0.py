import sympy

def generate_polynomial_string(n):
    """
    Generates the polynomial f_n(p) using the derived formula and returns it as a string.
    The formula is f_n(p) = (p^n - (1-p)^n) / (2p-1).
    This function uses sympy for symbolic manipulation.
    """
    if not isinstance(n, int) or n < 1:
        raise ValueError("n must be a positive integer.")

    p = sympy.Symbol('p')

    # The numerator of the formula
    numerator = p**n - (1 - p)**n
    # The denominator of the formula
    denominator = 2*p - 1
    
    # Sympy's expand function will perform the polynomial division and simplification
    poly = sympy.expand(numerator / denominator)
    
    return poly

def main():
    """
    Main function to demonstrate the formula by generating f_18(p) and printing its details.
    """
    # Let's generate f_18(p) as a verification example.
    n = 18
    f_n_p = generate_polynomial_string(n)

    # Pretty print the equation
    # We convert the sympy object to a string for printing
    poly_str = sympy.sstr(f_n_p, full_prec=False)
    
    print(f"The formula for f_n(p) is: (p**n - (1-p)**n) / (2*p - 1)\n")
    print(f"For n = {n}, the polynomial is:")
    print(f"f_{n}(p) = {poly_str}\n")
    
    # Per instruction, output each number in the final equation.
    # This is interpreted as printing the coefficients of the generated polynomial.
    # The coefficients are extracted from the highest power down to the constant term.
    coefficients = f_n_p.all_coeffs()
    
    print("The coefficients in the final equation, from the highest power of p to the constant term, are:")
    for coeff in coefficients:
        print(int(coeff))

if __name__ == '__main__':
    main()
