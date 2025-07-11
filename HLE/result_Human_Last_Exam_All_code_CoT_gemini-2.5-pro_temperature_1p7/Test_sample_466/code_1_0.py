# Plan:
# 1. Define the parameters n and q for GL(n, q).
# 2. Calculate each term in the product formula for the order of GL(n, q).
# 3. Compute the final product.
# 4. Print the calculation steps and the final answer, ensuring each number in the equation is shown.

# 1. Define parameters
n = 3
q = 11

# 2. Calculate the terms
term1 = q**n - q**0
term2 = q**n - q**1
term3 = q**n - q**2

# 3. Compute the product
order = term1 * term2 * term3

# 4. Print the detailed calculation and the final result
print(f"The number of internal adjunctions is the number of invertible 3x3 matrices over F_11, which is the order of the group GL(3, 11).")
print(f"The order is given by the formula: (q^n - q^0) * (q^n - q^1) * (q^n - q^2)")
print(f"For n={n} and q={q}, this evaluates to:")
print(f"({q}^{n} - {q**0}) * ({q}^{n} - {q**1}) * ({q}^{n} - {q**2})")
print(f"= ({q**n} - {q**0}) * ({q**n - q**1}) * ({q**n - q**2})")
print(f"= {term1} * {term2} * {term3}")
print(f"= {order}")