import math

def solve_poincare_polynomial():
    """
    Computes and prints the Poincaré polynomial for the given Lie algebra.
    """
    # Step 1: Define the dimensions of the Lie algebra and its commutator subalgebra.
    # The Lie algebra g is 6-dimensional.
    dim_g = 6

    # The commutator subalgebra [g, g] is spanned by the results of the non-zero brackets:
    # [e1, e2] = e3, [e1, e3] = e4, [e1, e5] = e6.
    # Thus, [g, g] = span{e3, e4, e6}.
    # Since e3, e4, e6 are basis vectors, they are linearly independent.
    dim_gg = 3

    # Step 2: Compute the dimension of the abelianization.
    # The Lie algebra is nilpotent. For a nilpotent Lie algebra, the homology is
    # determined by its abelianization g_ab = g/[g, g].
    # The dimension of the abelianization is dim(g) - dim([g, g]).
    dim_ab = dim_g - dim_gg

    # Step 3: Compute the coefficients of the Poincaré polynomial.
    # The Poincaré polynomial is P(x) = (1+x)^dim_ab.
    # The coefficient of x^k is the binomial coefficient C(dim_ab, k).
    coeffs = []
    for k in range(dim_ab + 1):
        coeffs.append(math.comb(dim_ab, k))

    # Step 4: Format and print the final polynomial, showing each number in the equation.
    # P(x) = c0*x^0 + c1*x^1 + c2*x^2 + ...
    poly_terms = []
    for k, c in enumerate(coeffs):
        # Format the term based on the power of x and the coefficient value.
        if c == 0:
            continue
        
        if k == 0:
            # Constant term
            term = str(c)
        elif k == 1:
            # Term with x
            if c == 1:
                term = "x"
            else:
                term = f"{c}*x"
        else:
            # Term with x^k
            if c == 1:
                term = f"x^{k}"
            else:
                term = f"{c}*x^{k}"
        poly_terms.append(term)

    # Join the terms with " + " to form the final polynomial string.
    poly_string = " + ".join(poly_terms)

    print("The Poincaré polynomial of the Lie algebra g is:")
    print(f"P(x) = {poly_string}")

solve_poincare_polynomial()