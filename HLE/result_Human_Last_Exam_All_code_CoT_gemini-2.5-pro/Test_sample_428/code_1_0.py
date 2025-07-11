# You can change these values to any positive integers
N1 = 1
M1 = 2
N2 = 3
M2 = 4

# The initial gaps between the particles are d1, d2, d3, d4
d1 = N1
d2 = M1
d3 = N2
d4 = M2

# The expected time is given by the formula:
# E[tau] = 1/3 * (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)

term1 = d1 * d2
term2 = d1 * d3
term3 = d1 * d4
term4 = d2 * d3
term5 = d2 * d4
term6 = d3 * d4

sum_of_products = term1 + term2 + term3 + term4 + term5 + term6
expected_tau = sum_of_products / 3

# The problem asks to output the final equation with numbers
# We will present the result in a readable format.
print(f"For N1={N1}, M1={M1}, N2={N2}, M2={M2}:")
print("The initial gaps are:")
print(f"d1 = {d1}, d2 = {d2}, d3 = {d3}, d4 = {d4}")
print("\nThe expectation of tau is E[tau] = 1/3 * (d1*d2 + d1*d3 + d1*d4 + d2*d3 + d2*d4 + d3*d4)")
print(f"E[tau] = 1/3 * ({d1}*{d2} + {d1}*{d3} + {d1}*{d4} + {d2}*{d3} + {d2}*{d4} + {d3}*{d4})")
print(f"E[tau] = 1/3 * ({term1} + {term2} + {term3} + {term4} + {term5} + {term6})")
print(f"E[tau] = 1/3 * ({sum_of_products})")
print(f"E[tau] = {expected_tau}")
