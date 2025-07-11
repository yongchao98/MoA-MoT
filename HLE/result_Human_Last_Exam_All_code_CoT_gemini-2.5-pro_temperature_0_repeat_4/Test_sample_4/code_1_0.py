import math

def compute_poincare_polynomial():
    """
    Computes and prints the Poincaré polynomial for the specified Lie algebra.

    The Lie algebra g is 6-dimensional with relations:
    [e1, e2] = e3, [e1, e3] = e4, [e1, e5] = e6.

    This algebra is nilpotent. For a nilpotent Lie algebra, the Poincaré
    polynomial is determined by the dimension of its abelianization, b1.
    The Betti numbers are b_k = C(b1, k).
    """
    # Dimension of the Lie algebra g
    dim_g = 6

    # The derived algebra [g, g] is spanned by the results of the commutators.
    # The basis of [g, g] is {e3, e4, e6}.
    dim_derived_g = 3

    # The dimension of the first homology group H_1(g), denoted as b1, is
    # the dimension of the abelianization g / [g, g].
    b1 = dim_g - dim_derived_g

    # The Betti numbers b_k = dim(H_k(g)) are given by the binomial
    # coefficients C(b1, k).
    betti_numbers = []
    for k in range(b1 + 1):
        betti_numbers.append(math.comb(b1, k))

    # Construct the Poincaré polynomial string P(x) = sum(b_k * x^k).
    # We ensure every coefficient is explicitly printed.
    poly_terms = []
    for k, coeff in enumerate(betti_numbers):
        if k == 0:
            term = str(coeff)
        elif k == 1:
            term = f"{coeff}*x"
        else:
            term = f"{coeff}*x^{k}"
        poly_terms.append(term)
    
    polynomial_string = " + ".join(poly_terms)
    
    print(f"The Poincaré polynomial is P(x) = {polynomial_string}")

compute_poincare_polynomial()