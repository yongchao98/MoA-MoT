# Set the parameters for the finite field and vector space dimension
q = 11
n = 3

# Calculate the powers of q needed for the formula
q_n = q**n
q_0 = q**0
q_1 = q**1
q_2 = q**2

# Calculate each term in the product for |GL(n, F_q)|
term1 = q_n - q_0
term2 = q_n - q_1
term3 = q_n - q_2

# Calculate the final result
result = term1 * term2 * term3

# Print the step-by-step calculation
print(f"The number of adjunctions is the size of the general linear group GL(n, F_q).")
print(f"The formula is: |GL(n, F_q)| = (q^n - q^0) * (q^n - q^1) * ... * (q^n - q^(n-1))")
print(f"For this problem, n = {n} and q = {q}.")
print(f"The equation is: ({q}^{n} - {q}^{0}) * ({q}^{n} - {q}^{1}) * ({q}^{n} - {q}^{2})")
print(f"= ({q_n} - {q_0}) * ({q_n} - {q_1}) * ({q_n} - {q_2})")
print(f"= {term1} * {term2} * {term3}")
print(f"= {result}")
