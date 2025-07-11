import sympy

def solve_and_sum_coordinates():
    """
    This function solves the given matrix equations to find the first coordinates
    of the solution matrices and calculates their sum.
    """
    # Define symbols for symbolic computation
    i = sympy.I
    sqrt = sympy.sqrt

    # --- Equation 1 ---
    # From the analysis, we have a^2 = -53/132 for the first coordinate of X1.
    # The solutions for the first coordinate of X1 are:
    val1_sq = sympy.Rational(-53, 132)
    a1_1 = sqrt(val1_sq)
    a1_2 = -sqrt(val1_sq)

    # --- Equation 2 ---
    # From the analysis, we have a^2 = -3/110 for the first coordinate of X2.
    # The solutions for the first coordinate of X2 are:
    val2_sq = sympy.Rational(-3, 110)
    a2_1 = sqrt(val2_sq)
    a2_2 = -sqrt(val2_sq)

    # The problem asks for the sum of the first coordinates of all solutions.
    total_sum = a1_1 + a1_2 + a2_1 + a2_2

    # Print the final equation showing each term and the result.
    # The terms are formatted to show the symbolic representation.
    term1_str = f"({sympy.printing.sstr(a1_1)})"
    term2_str = f"({sympy.printing.sstr(a1_2)})"
    term3_str = f"({sympy.printing.sstr(a2_1)})"
    term4_str = f"({sympy.printing.sstr(a2_2)})"
    
    print("The sum of the first coordinate of solutions is calculated as follows:")
    print(f"{term1_str} + {term2_str} + {term3_str} + {term4_str} = {total_sum}")

if __name__ == "__main__":
    solve_and_sum_coordinates()