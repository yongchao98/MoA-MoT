import sympy as sp

def main():
    """
    This function computes the total Chern class of the bundle S^(2,1)T(CP^4).
    The plan is as follows:
    1. Let E = T(CP^4). Its total Chern class is c(E) = (1+h)^5, where h is the hyperplane class.
       This means the elementary symmetric polynomials of the formal Chern roots x_i of E are
       e_k(x) = C(5, k) * h^k.
    2. We use the splitting principle. The Chern roots of S^(2,1)E can be determined from the
       character of the representation S^(2,1) of GL(4). The weights are of the form
       2*x_i + x_j (for i != j) and x_i + x_j + x_k (for distinct i,j,k, with multiplicity 2).
    3. The total Chern class of S^(2,1)E is c = 1 + c_1 + c_2 + ...
       The Chern classes c_k are the elementary symmetric polynomials of these weights.
    4. We can express c_k in terms of the power sums P_j = sum(w_i^j) of the weights, using Newton's identities.
    5. Each P_j(w) is a symmetric polynomial in the original roots x_i. We express them
       in the basis of power sums P_k(x) = sum(x_i^k).
    6. Finally, we use Newton's identities again to relate P_k(x) back to e_k(x) = C(5, k)*h^k.
    7. All computations are done modulo h^5, as h^5=0 in the cohomology ring of CP^4.
    """
    h = sp.Symbol('h')
    N = 4 # for CP^4

    # Formal roots of T(CP^4)
    roots_E = sp.symbols(f'x_1:{N + 1}')

    # Power sums for these roots, as expressions in h
    p_x_h = [sp.Integer(0)] * (N + 1)
    e_x_h = [sp.binomial(N + 1, k) * h**k for k in range(N + 1)]
    p_x_h[1] = e_x_h[1]
    for k in range(2, N + 1):
        p_x_h[k] = (-1)**(k-1) * k * e_x_h[k] + sum([(-1)**(j-1) * e_x_h[j] * p_x_h[k-j] for j in range(1, k)])

    def truncate_h(expr):
        """Truncate expression at h^5"""
        return sp.poly(expr, h).trunc(N + 1).as_expr()

    for k in range(1, N + 1):
        p_x_h[k] = truncate_h(p_x_h[k])

    # Weights of S^(2,1)E
    weights = []
    from itertools import combinations, permutations
    for i, j in permutations(range(N), 2):
        weights.append(2 * roots_E[i] + roots_E[j])
    for i, j, k in combinations(range(N), 3):
        w = roots_E[i] + roots_E[j] + roots_E[k]
        weights.append(w)
        weights.append(w)

    # Power sums of roots_E, as polynomials
    p_x_poly = [sp.Integer(0)]*(N+1)
    for k in range(1, N+1):
        p_x_poly[k] = sum(r**k for r in roots_E)

    # Power sums of the weights, expressed in terms of p_x_poly
    p_w_poly = [sp.Integer(0)] * (N + 1)
    for k in range(1, N + 1):
        pk_w_raw = sum(w**k for w in weights)

        from sympy.polys.polyfuncs import symmetrize
        # symmetrize returns a tuple (symmetric_poly, G), we take the poly
        # The result is in terms of elementary symmetric polynomials e_k
        e_syms = sp.symbols('e1:{}'.format(N+1))
        # symmetrize does not work for higher degree with this many terms.
        # Fallback to linear algebra
        from sympy.utilities.iterables import partitions
        basis_parts = partitions(k)
        basis = []
        for part in basis_parts:
            term = 1
            for p_idx, p_pow in part.items():
                term *= p_x_poly[p_idx]**p_pow
            basis.append(term)
        
        coeffs = sp.symbols(f'a_1:{len(basis)+1}')
        generic_form = sum(c * b for c, b in zip(coeffs, basis))
        
        eqs = []
        num_trials = len(coeffs)
        for i in range(num_trials):
            rand_vals = [i + l + 2 for l in range(N)] # Use deterministic values to ensure reproducibility
            subs_map = {r: v for r, v in zip(roots_E, rand_vals)}
            raw_val = pk_w_raw.subs(subs_map)
            form_val = generic_form.subs(subs_map)
            eqs.append(sp.Eq(form_val, raw_val))

        sol = sp.solve(eqs, coeffs, dict=True)
        if sol:
             p_w_poly[k] = generic_form.subs(sol[0])
        else:
             print(f"Error: could not solve for P_{k}(w)")
             p_w_poly[k] = 0

    # Chern classes of S^(2,1)E in terms of power sums of weights (Newton's identities)
    c_w_poly = [sp.Integer(0)] * (N + 1)
    c_w_poly[0] = sp.Integer(1)
    c_w_poly[1] = p_w_poly[1]
    for k in range(2, N + 1):
        c_w_poly[k] = ((-1)**(k-1) * p_w_poly[k] + sum([(-1)**(j-1) * c_w_poly[k-j] * p_w_poly[j] for j in range(1, k)])) / k
        
    # Final step: substitute p_k(x) expressions in h
    subs_map = {p_x_poly[k]: p_x_h[k] for k in range(1, N + 1)}
    
    final_c = [sp.Integer(0)] * (N + 1)
    final_c[0] = sp.Integer(1)
    for k in range(1, N + 1):
        final_c[k] = truncate_h(sp.expand(c_w_poly[k].subs(subs_map)))

    total_chern = sum(final_c)
    
    # The problem asks for the expression as a polynomial in h
    # and to output the equation in the final response.
    # We will build the string representation of the final answer.
    
    pretty_total_chern = sp.poly(total_chern, h).as_expr()
    final_coeffs = sp.poly(pretty_total_chern, h).all_coeffs()
    if len(final_coeffs) < N+1:
        final_coeffs = [0]*(N+1-len(final_coeffs)) + final_coeffs
    final_coeffs.reverse() # from c0 to c4

    c1, c2, c3, c4 = final_coeffs[1:]
    
    print(f"c(S^(2,1)} T(CP^4)) = 1 + ({c1}) h + ({c2}) h^2 + ({c3}) h^3 + ({c4}) h^4")
    
if __name__ == '__main__':
    main()