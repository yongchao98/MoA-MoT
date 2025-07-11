import spherogram

def get_knot_jones_polynomial():
    """
    This function calculates the Jones polynomial for the knot shown in the image (10_124)
    and prints it in a standard mathematical format.
    """
    try:
        # The knot in the image is the 10_124 knot, also known as a Perko pair knot.
        # We can create a Knot object for it using its standard name.
        knot = spherogram.Knot('10_124')

        # The spherogram library calculates a specific version of the Jones polynomial.
        # Let's call it J(t). It is usually a polynomial, not a Laurent polynomial.
        # The variable used by the library is symbolic, typically 'q' or 't'.
        jones_poly_J = knot.jones_polynomial()

        # The desired format is a Laurent polynomial V(t), which is related to J(t)
        # by the change of variable t -> t^-1. So, V(t) = J(t^-1).
        # We can get the polynomial as a dictionary of {power: coefficient}.
        poly_dict_J = jones_poly_J.dict()

        # To get V(t), we negate the powers of J(t).
        poly_dict_V = {-power: coeff for power, coeff in poly_dict_J.items()}

        # Sort the powers in descending order for formatting.
        sorted_powers = sorted(poly_dict_V.keys(), reverse=True)

        # Build the polynomial string term by term.
        result_terms = []
        for i, power in enumerate(sorted_powers):
            coeff = poly_dict_V[power]
            if coeff == 0:
                continue

            # Determine the sign prefix for the term.
            if i == 0:
                sign = "-" if coeff < 0 else ""
            else:
                sign = " - " if coeff < 0 else " + "

            abs_coeff = abs(coeff)

            # Format the coefficient part of the term.
            # By convention, a coefficient of 1 is omitted unless it's a constant term.
            if abs_coeff == 1 and power != 0:
                coeff_str = ""
            else:
                coeff_str = str(abs_coeff)

            # Format the variable and power part of the term.
            if power == 0:
                var_str = ""
            elif power == 1:
                var_str = "t"
            else:
                # Use curly braces for multi-character exponents for clarity.
                var_str = f"t^{{{power}}}"
            
            result_terms.append(f"{sign}{coeff_str}{var_str}")

        # Join all terms and print the final polynomial string.
        final_polynomial = "".join(result_terms)
        print(final_polynomial)

    except ImportError:
        print("The 'spherogram' library is required. Please install it using 'pip install spherogram'.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    get_knot_jones_polynomial()