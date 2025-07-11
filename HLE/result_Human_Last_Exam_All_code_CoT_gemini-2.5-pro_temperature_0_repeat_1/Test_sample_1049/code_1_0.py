import math

def nCr_safe(n, r):
    """
    Calculates binomial coefficient C(n, r), returns 0 if n < r.
    """
    if r < 0 or r > n:
        return 0
    r = min(r, n - r)
    if r == 0:
        return 1
    numer = 1
    denom = 1
    for i in range(r):
        numer *= (n - i)
        denom *= (i + 1)
    return numer // denom

def get_poly_coeffs():
    """
    Calculates the coefficients of the polynomial Q(n) in the binomial basis.
    Q(n) is defined by S_n = 4^n * Q(n).
    The coefficients c_k are determined by c_k = Delta^k Q(0), where Delta is the forward difference operator.
    """
    p = [1, 464, 8256, 18944, 4096]
    
    q_values = []
    for n in range(6):
        s_n = 0
        for j in range(min(n, 4) + 1):
            s_n += p[j] * nCr_safe(n - j + 5, 5) * (4**(n - j))
        q_values.append(s_n // (4**n))

    # Calculate forward differences at n=0 to find coefficients
    coeffs = []
    diff_table = [list(q_values)]
    for k in range(5):
        prev_diff = diff_table[-1]
        new_diff = [prev_diff[i+1] - prev_diff[i] for i in range(len(prev_diff)-1)]
        diff_table.append(new_diff)

    coeffs.append(diff_table[0][0]) # c0 = Q(0)
    coeffs.append(diff_table[1][0]) # c1 = Delta Q(0)
    coeffs.append(diff_table[2][0]) # c2 = Delta^2 Q(0)
    coeffs.append(diff_table[3][0]) # c3 = Delta^3 Q(0)
    coeffs.append(diff_table[4][0]) # c4 = Delta^4 Q(0)
    coeffs.append(diff_table[5][0]) # c5 = Delta^5 Q(0)
    
    return coeffs

def main():
    """
    Prints the closed form for the sum.
    """
    coeffs = get_poly_coeffs()
    
    # Reverse coeffs for printing from highest degree
    coeffs.reverse()
    
    formula_parts = []
    for i, c in enumerate(coeffs):
        degree = 5 - i
        if c != 0:
            if degree > 1:
                formula_parts.append(f"{c} * C(n, {degree})")
            elif degree == 1:
                formula_parts.append(f"{c} * n")
            else: # degree == 0
                formula_parts.append(f"{c}")

    final_formula = "4^n * (" + " + ".join(formula_parts) + ")"
    
    print("The closed form for the sum is:")
    print(final_formula)

if __name__ == "__main__":
    main()
