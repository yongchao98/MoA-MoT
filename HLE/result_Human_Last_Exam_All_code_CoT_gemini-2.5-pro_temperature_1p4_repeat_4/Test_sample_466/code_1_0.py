q = 11
n = 3

q_n = q**n

term1 = q_n - 1
term2 = q_n - q
term3 = q_n - q**2

result = term1 * term2 * term3

print(f"The number of internal adjunctions is the order of GL_3(F_11), which is calculated as:")
print(f"({q**n} - {q**0}) * ({q**n} - {q**1}) * ({q**n} - {q**2})")
print(f"= {term1} * {term2} * {term3}")
print(f"= {result}")