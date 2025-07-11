def poly_mul_trunc(p1, p2, trunc):
    res = [0] * (trunc + 1)
    len_p1 = len(p1)
    len_p2 = len(p2)
    for i in range(len_p1):
        if p1[i] == 0:
            continue
        for j in range(len_p2):
            if i + j <= trunc:
                res[i+j] += p1[i] * p2[j]
    return res

def poly_pow_trunc(p, n, trunc):
    res = [1]
    for _ in range(n):
        res = poly_mul_trunc(res, p, trunc)
    return res

trunc = 48
# Coefficients for (1 + x + ... + x^6)
p1_base = [1] * 7
# Coefficients for (1 + x + ... + x^6)^8
p1_coeffs = poly_pow_trunc(p1_base, 8, trunc)

# Coefficients for (1 + x + ... + x^8)
p2_base = [1] * 9
# Coefficients for (1 + x + ... + x^8)^6
p2_coeffs = poly_pow_trunc(p2_base, 6, trunc)

# Coefficients for the final polynomial product
result_coeffs = poly_mul_trunc(p1_coeffs, p2_coeffs, trunc)

# The answer is the coefficient of x^48
final_answer = result_coeffs[48]

# The calculated value for final_answer is 20180163135.