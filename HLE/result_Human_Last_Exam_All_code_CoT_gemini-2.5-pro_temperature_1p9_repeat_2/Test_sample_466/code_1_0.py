import math

q = 11
n = 3

# Calculate the terms of the product formula for the order of GL(n, q)
term1 = q**n - 1
term2 = q**n - q
term3 = q**n - q**2

# Calculate the total number of adjunctions
num_adjunctions = term1 * term2 * term3

# Print the calculation step-by-step
print(f"The number of internal adjunctions is the order of the general linear group GL({n}, F_{q}).")
print(f"The formula is |GL({n}, F_{q})| = (q^n - q^0) * (q^n - q^1) * ... * (q^n - q^(n-1))")
print(f"For q={q} and n={n}:")
print(f"Number of adjunctions = ({q}^{n} - {q**0}) * ({q}^{n} - {q**1}) * ({q}^{n} - {q**2})")
print(f"Number of adjunctions = ({q**n} - {1}) * ({q**n} - {q}) * ({q**n} - {q**2})")
print(f"Number of adjunctions = ({term1}) * ({term2}) * ({term3})")
print(f"Final calculation: {term1} * {term2} * {term3} = {num_adjunctions}")
print(f"\nThe total number of internal adjunctions is: {num_adjunctions}")
