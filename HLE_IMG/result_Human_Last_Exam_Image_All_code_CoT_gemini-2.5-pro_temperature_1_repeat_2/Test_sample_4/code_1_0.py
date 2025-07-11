import sys
import subprocess

# This script requires the 'spherogram' and 'sympy' libraries.
# The following code will attempt to install them if they are not found.
try:
    import spherogram
    import sympy
except ImportError:
    print("Required libraries not found. Attempting to install 'spherogram' and 'sympy'...")
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "spherogram", "sympy"])
        import spherogram
        import sympy
        print("Libraries installed successfully.")
    except Exception as e:
        print(f"Error: Failed to install required libraries. Please install them manually using 'pip install spherogram sympy'.")
        print(f"Details: {e}")
        sys.exit(1)

def get_knot_jones_polynomial_str(pd_code):
    """
    Calculates the Jones polynomial for a knot given its PD code and formats it as a string.

    Args:
        pd_code (list): The Planar Diagram (PD) code of the knot.

    Returns:
        str: The formatted Jones polynomial string.
    """
    # Create a Link object from the PD code
    link = spherogram.Link(pd_code)
    
    # Define the symbolic variable for the polynomial
    t = sympy.Symbol('t')
    
    # Calculate the Jones polynomial
    jones_poly_expr = link.jones_polynomial(t)

    # --- Format the polynomial into the desired string format ---

    # Extract the terms and their coefficients from the sympy expression
    terms_dict = jones_poly_expr.as_coefficients_dict()

    # Create a list of (degree, coefficient) tuples
    term_list = []
    for term_key, coeff in terms_dict.items():
        if term_key.is_number:
            degree = 0
        else:
            degree = sympy.degree(term_key, t)
        
        # Handle terms like 1/t^n which sympy might represent with degree 1
        if term_key.is_Pow and term_key.exp < 0:
             degree = term_key.exp

        term_list.append((int(degree), int(coeff)))

    # Sort terms by degree in descending order
    term_list.sort(key=lambda x: x[0], reverse=True)

    # Build the output string part by part
    poly_str_parts = []
    for i, (degree, coeff) in enumerate(term_list):
        # 1. Determine the sign and spacing
        sign_str = ""
        if i == 0:
            if coeff < 0:
                sign_str = "-"
        else:
            if coeff > 0:
                sign_str = " + "
            else:
                sign_str = " - "
        
        # 2. Determine the coefficient string
        coeff_val = abs(coeff)
        coeff_str = ""
        # Only show coefficient if it's not 1, or if it's the constant term
        if coeff_val != 1 or degree == 0:
            coeff_str = str(coeff_val)
            
        # 3. Determine the variable part string
        var_str = ""
        if degree != 0:
            # Add 't' if the coefficient is not shown (i.e., it's 1 or -1)
            if coeff_val == 1:
                var_str = "t"
            # Add '*t' if coefficient is shown
            else:
                var_str = "t" # Style is 5t, not 5*t
        
        if degree != 1 and degree != 0:
            # Use {exponent} for clarity, as in the example t^{-10}
            var_str += f"^{{{degree}}}"
            
        # 4. Combine parts for the current term
        poly_str_parts.append(f"{sign_str}{coeff_str}{var_str}")

    return "".join(poly_str_parts)

if __name__ == "__main__":
    # The knot in the image is the 8_19 knot.
    # Its Planar Diagram (PD) code can be found in knot theory resources like Knot Atlas.
    pd_code_8_19 = [
        [1, 5, 2, 6], [7, 11, 8, 12], [9, 1, 10, 2], [3, 9, 4, 10],
        [13, 3, 14, 4], [5, 15, 6, 16], [11, 13, 12, 14], [15, 7, 16, 8]
    ]
    
    # Get the formatted Jones polynomial
    jones_polynomial_string = get_knot_jones_polynomial_str(pd_code_8_19)
    
    # Print the final result
    print("The Jones polynomial p(t) of the knot is:")
    print(jones_polynomial_string)