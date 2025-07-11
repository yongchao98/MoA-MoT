import mpmath

def solve_machin_like_formula():
    """
    Finds the integer coefficients for a given Machin-like formula for pi.
    This script requires the 'mpmath' library. You can install it using:
    pip install mpmath
    """
    try:
        # Set high precision for mpmath for the PSLQ algorithm to work reliably.
        mpmath.mp.dps = 200

        # These are the denominators from the problem statement.
        denominators = [122, 239, 682, 1252, 2855, 12943]

        # We are looking for an integer relation for the following set of real numbers:
        # arctan(1/122), arctan(1/239), ..., pi/4
        terms = [mpmath.atan(mpmath.mpf(1)/d) for d in denominators]
        terms.append(mpmath.pi / 4)

        # Use the PSLQ algorithm to find an integer relation.
        # It finds a non-zero integer vector 'coeffs' such that sum(coeffs[i] * terms[i]) is close to 0.
        coeffs = mpmath.pslq(terms)

        # The relation found is of the form:
        # c_1*atan(1/x_1) + ... + c_6*atan(1/x_6) + k*(pi/4) = 0
        # The problem asks for the form:
        # c_1*atan(1/x_1) + ... = n*(pi/4)
        # Comparing these gives n = -k, where k is the last coefficient found by PSLQ.
        k = coeffs[-1]

        # The problem requires the smallest *positive* n.
        # PSLQ finds a primitive relation, so abs(n) is the smallest possible.
        # If n is negative (i.e., k is positive), we multiply the entire relation by -1.
        if k > 0:
            coeffs = [-c for c in coeffs]
        
        # Extract the coefficients c_1, ..., c_6 and n.
        c = coeffs[:-1]
        n = -coeffs[-1]

        # As requested, output each number in the final equation.
        print("The discovered relation is:")
        
        equation_lhs = f"{n} * pi / 4"
        
        equation_rhs_parts = []
        is_first_term = True
        for i in range(len(c)):
            if c[i] != 0:
                sign = "+" if c[i] > 0 else "-"
                abs_coeff = abs(c[i])
                
                term = f"{abs_coeff}*arctan(1/{denominators[i]})"
                
                if is_first_term:
                    if sign == "-":
                        equation_rhs_parts.append(f"-{term}")
                    else:
                        equation_rhs_parts.append(term)
                    is_first_term = False
                else:
                    equation_rhs_parts.append(f" {sign} {term}")

        equation_rhs = "".join(equation_rhs_parts)
        print(f"{equation_lhs} = {equation_rhs}")
        print("\n")

        # Provide the answer in the specified format for the final result.
        result_string = f"{n},{','.join(map(str, c))}"
        print("The solution in the format n,c1,c2,c3,c4,c5,c6 is:")
        print(result_string)

    except ImportError:
        print("This script requires the 'mpmath' library.")
        print("Please install it using: pip install mpmath")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == '__main__':
    solve_machin_like_formula()