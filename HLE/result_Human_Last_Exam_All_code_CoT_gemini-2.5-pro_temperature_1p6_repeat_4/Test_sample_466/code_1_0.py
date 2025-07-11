q = 11
n = 3

q_pow_n = q**n
q_pow_0 = q**0
q_pow_1 = q**1
q_pow_2 = q**2

term1 = q_pow_n - q_pow_0
term2 = q_pow_n - q_pow_1
term3 = q_pow_n - q_pow_2

result = term1 * term2 * term3

print(f"The number of internal adjunctions is the order of GL_3(F_11), calculated as:")
print(f"({q}^3 - {q}^0) * ({q}^3 - {q}^1) * ({q}^3 - {q}^2) = ({q_pow_n} - {q_pow_0}) * ({q_pow_n} - {q_pow_1}) * ({q_pow_n} - {q_pow_2})")
print(f"= {term1} * {term2} * {term3}")
print(f"= {result}")