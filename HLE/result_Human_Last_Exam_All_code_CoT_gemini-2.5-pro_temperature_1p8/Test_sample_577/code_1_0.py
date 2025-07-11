from sympy import symbols, Poly, series

def chern_to_power_sums(c, n, h, ring):
    """Converts Chern classes to power sums."""
    p = [ring(0)] * (n + 1)
    if n >= 1:
        p[1] = c[1]
    for k in range(2, n + 1):
        s = c[1] * p[k-1]
        for i in range(2, k):
            s -= c[i] * p[k-i]
        s += k * c[k]
        p[k] = s
    return p

def power_sums_to_chern(p, n, h, ring):
    """Converts power sums to Chern classes."""
    c = [ring(1)] + [ring(0)] * n
    if n >= 1:
        c[1] = p[1]
    for k in range(2, n + 1):
        s = ring(0)
        for i in range(1, k):
            s -= c[i] * p[k-i]
        s -= p[k]
        c[k] = s / (-k)
    return c
    
def main():
    h = symbols('h')
    
    # Cohomology ring H*(CP^4) is Z[h]/(h^5)
    # We use sympy's series to truncate polynomials at h^5
    def make_poly(coeffs):
        p = 0
        for i, c in enumerate(coeffs):
            p += c * h**i
        return series(p, h, 0, 5).as_poly(h)

    def ring(p):
        return series(p, h, 0, 5).as_poly(h)

    # Chern classes of E = TCP^4
    # c(E) = (1+h)^5 = 1 + 5h + 10h^2 + 10h^3 + 5h^4
    c_E = [ring(1), ring(5*h), ring(10*h**2), ring(10*h**3), ring(5*h**4)]
    rank_E = 4

    # Power sums of E
    p_E = chern_to_power_sums(c_E, rank_E, h, ring)
    
    # Analysis of Lambda^2(E)
    rank_L2E = 6
    p_L2E = [ring(0)] * (rank_E + 1)
    p_L2E[1] = (rank_E - 1) * p_E[1]
    p_L2E[2] = (rank_E - 2) * p_E[2] + p_E[1]**2 - p_E[2]  # (n-2)p_2 + 2e_2
    p_L2E[3] = (rank_E - 3) * p_E[3] + 3*(p_E[1]*p_E[2] - p_E[3]) # (n-3)p_3 + 3e_2p_1 - 3e_3, or simpler forms.
    p_L2E[4] = (rank_E - 4) * p_E[4] + p_E[1]*p_E[3] - p_E[4] + p_E[2]**2 - p_E[4] + 3*(p_E[1]**2*p_E[2] - p_E[1]**2*p_E[2]...)
    
    # simpler rederived formulas for p_k(Lambda^2 E)
    p_L2E[2] = 3 * p_E[2] + 2 * c_E[2] # 3*p_2 + 2*e_2
    p_L2E[3] = 3 * p_E[1] * p_E[2] # n=4 specific
    p_L2E[4] = 3 * p_E[4] + 4*p_E[3]*p_E[1] + p_E[2]**2 # n=4 specific
    c_L2E = power_sums_to_chern(p_L2E, rank_L2E, h, ring)

    # Analysis of Lambda^3(E)
    rank_L3E = 4
    p_L3E = [ring(0)] * (rank_E + 1)
    p_L3E[1] = Poly(3*5*h, h) # 3*c1(E)
    y_roots = [c_E[1] - r for r in [0,0,0,0]] # formal roots for calculation
    # simplified expressions for p_k(Lambda^3 E) from roots e_1 - x_i
    p_L3E[1] = 4*c_E[1] - p_E[1]
    p_L3E[2] = 4*c_E[1]**2 - 2*c_E[1]*p_E[1] + p_E[2]
    p_L3E[3] = 4*c_E[1]**3 - 3*c_E[1]**2*p_E[1] + 3*c_E[1]*p_E[2] - p_E[3]
    p_L3E[4] = 4*c_E[1]**4 - 4*c_E[1]**3*p_E[1] + 6*c_E[1]**2*p_E[2] - 4*c_E[1]*p_E[3] + p_E[4]
    c_L3E = power_sums_to_chern(p_L3E, rank_L3E, h, ring)
    
    # Analysis of E tensor Lambda^2(E)
    F = E
    G = {'rank': rank_L2E, 'p': p_L2E, 'c': c_L2E}
    
    p_FG = [ring(0)] * 5
    p_FG[1] = F['rank']*G['p'][1] + G['rank']*F['p'][1]
    p_FG[2] = F['rank']*G['p'][2] + G['rank']*F['p'][2] + 2*F['p'][1]*G['p'][1]
    p_FG[3] = F['rank']*G['p'][3] + G['rank']*F['p'][3] + 3*(F['p'][1]*G['p'][2] + F['p'][2]*G['p'][1])
    p_FG[4] = F['rank']*G['p'][4] + G['rank']*F['p'][4] + 4*(F['p'][1]*G['p'][3] + F['p'][3]*G['p'][1]) + 6*F['p'][2]*G['p'][2]
    c_FG = power_sums_to_chern(p_FG, 4, h, ring)

    # The actual bundles
    F = {'name': 'E', 'rank': rank_E, 'p': p_E, 'c': c_E}
    G = {'name': 'L2E', 'rank': rank_L2E, 'p': p_L2E, 'c': c_L2E}
    
    p_FG = [ring(0)] * 5
    p_FG[1] = F['rank']*G['p'][1] + G['rank']*F['p'][1]
    p_FG[2] = F['rank']*G['p'][2] + G['rank']*F['p'][2] + 2*F['p'][1]*G['p'][1]
    p_FG[3] = F['rank']*G['p'][3] + G['rank']*F['p'][3] + 3*(F['p'][1]*G['p'][2] + F['p'][2]*G['p'][1])
    p_FG[4] = F['rank']*G['p'][4] + G['rank']*F['p'][4] + 4*(F['p'][1]*G['p'][3] + F['p'][3]*G['p'][1]) + 6*F['p'][2]*G['p'][2]
    
    # We only need Chern classes up to degree 4. Rank of tensor product is 24.
    c_tensor = power_sums_to_chern(p_FG, 4, h, ring)

    # Compute c(S^{(2,1)}E) = c(E tensor L2E) / c(L3E)
    C_tensor = Poly(reversed(c_tensor), h)
    C_L3E = Poly(reversed(c_L3E), h)
    
    # Polynomial division
    Q, R = C_tensor.div(C_L3E)

    # The result needs to be printed in the correct order for a Chern class polynomial.
    coeffs = Q.all_coeffs()
    coeffs.reverse()
    
    final_c = "1"
    for i in range(1, 5):
      if i < len(coeffs):
        term_coeff = coeffs[i].coeff(h,i)
        if term_coeff != 0:
            final_c += f" + {term_coeff}*h^{i}"

    print(f"c(S^(2,1) T(CP^4)) = 1 + {coeffs[1].coeff(h,1)}*h + {coeffs[2].coeff(h,2)}*h^2 + {coeffs[3].coeff(h,3)}*h^3 + {coeffs[4].coeff(h,4)}*h^4")

main()