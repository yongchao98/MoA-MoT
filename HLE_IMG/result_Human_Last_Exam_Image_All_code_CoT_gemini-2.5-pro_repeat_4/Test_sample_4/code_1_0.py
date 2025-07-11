import spherogram

def get_jones_polynomial_for_10_139():
    """
    Calculates and prints the Jones polynomial for the knot 10_139.

    The knot in the image is identified as 10_139. This function uses the
    spherogram library to compute its Jones polynomial. The knot is defined
    by its Planar Diagram (PD) code. The resulting polynomial is then
    printed in decreasing degree order, with each component of the
    equation printed separately as requested.
    """
    # The Planar Diagram (PD) code for the knot 10_139.
    # This is a standard representation of the knot's projection.
    pd_code = [[8, 10, 9, 11], [18, 2, 19, 1], [20, 16, 1, 15], [16, 4, 17, 3],
               [4, 14, 5, 13], [14, 12, 15, 11], [12, 6, 13, 7], [6, 18, 7, 17],
               [2, 8, 3, 9], [10, 20, 11, 19]]

    # Create a Link object from the PD code.
    knot = spherogram.Link(pd_code)

    # Calculate the Jones polynomial. The result is a LaurentPolynomial object.
    p = knot.jones_polynomial()

    # Get the terms of the polynomial as a dictionary of {degree: coefficient}.
    poly_terms = p.dict

    # Sort the terms by degree in descending order.
    sorted_terms = sorted(poly_terms.items(), key=lambda item: item[0], reverse=True)

    is_first_term = True
    for degree, coeff in sorted_terms:
        # Determine the sign and print it.
        sign = ""
        if not is_first_term:
            if coeff > 0:
                sign = " + "
            else:
                sign = " - "
        elif coeff < 0:
            sign = "-"
        
        print(sign, end="")

        # Get the absolute value of the coefficient for printing.
        abs_coeff = abs(coeff)

        # Print the coefficient, unless it's 1 and not a constant term.
        if abs_coeff != 1 or degree == 0:
            print(abs_coeff, end="")

        # Print the variable 't' and its exponent, if applicable.
        if degree != 0:
            print("t", end="")
            if degree != 1:
                print(f"^{degree}", end="")
        
        is_first_term = False
    
    print() # For a final newline

get_jones_polynomial_for_10_139()
