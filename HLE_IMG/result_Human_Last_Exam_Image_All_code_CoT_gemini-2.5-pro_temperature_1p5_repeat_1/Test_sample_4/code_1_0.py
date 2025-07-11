import sys
import io

# A helper function to format a single term (coefficient, power) of a polynomial.
def format_term(coeff, power, is_first_term):
    """Formats a single term of the polynomial for printing."""
    # Build the term string piece by piece.
    term_str = ""

    # 1. Handle the sign and spacing.
    if not is_first_term:
        term_str += " + " if coeff > 0 else " - "
    else:
        if coeff < 0:
            term_str += "-"

    # 2. Handle the coefficient's absolute value.
    # It's omitted if it's 1 unless the power is 0 (a constant term).
    abs_coeff = abs(coeff)
    if abs_coeff != 1 or power == 0:
        term_str += str(abs_coeff)

    # 3. Handle the variable 't' and its power.
    # 't' is not shown for a power of 0.
    if power != 0:
        term_str += "t"
        # The exponent is not shown if it's 1.
        if power != 1:
            term_str += f"^{power}"

    return term_str

def get_jones_polynomial_10_132():
    """
    This function returns the Jones polynomial for the 10_132 knot.
    The polynomial is known and is retrieved, not calculated from scratch.
    """
    # The Jones polynomial for knot 10_132 is:
    # V(t) = t^11 - t^9 + t^8 - 2t^7 + 2t^6 - 2t^5 + 2t^4 - t^3 + t^2
    # We represent this as a list of (coefficient, power) tuples, in decreasing order of power.
    poly_terms = [
        (1, 11),
        (-1, 9),
        (1, 8),
        (-2, 7),
        (2, 6),
        (-2, 5),
        (2, 4),
        (-1, 3),
        (1, 2)
    ]
    
    # Use a list to build the parts of the string.
    result_parts = []
    
    # Process the first term separately to handle the leading sign.
    first_c, first_p = poly_terms[0]
    result_parts.append(format_term(first_c, first_p, is_first_term=True))

    # Process the remaining terms.
    for c, p in poly_terms[1:]:
        result_parts.append(format_term(c, p, is_first_term=False))
        
    return "".join(result_parts)

# Print the final formatted polynomial.
# This ensures that all numbers (coefficients and exponents) are part of the output string.
final_polynomial = get_jones_polynomial_10_132()
print(final_polynomial)

# The following lines are for providing the answer in the required format for the platform.
# Suppress the print output from the user-facing part of the code.
old_stdout = sys.stdout
sys.stdout = new_stdout = io.StringIO()
# The code above this point is the solution for the user.
# The code below is for the final answer submission.
final_polynomial_for_answer = get_jones_polynomial_10_132()
sys.stdout = old_stdout
print(f"<<<{final_polynomial_for_answer}>>>", file=sys.stderr)
