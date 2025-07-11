#
# Step 1: Define the parameters of our problem.
# We are working with the vector space F_11^3.
# So n=3 and q=11.
n = 3
q = 11

# Step 2: Calculate the terms needed for the GL(n, q) formula.
# The terms are (q^n - q^0), (q^n - q^1), ..., (q^n - q^(n-1)).

# Calculate q^n
q_n = q**n

# The terms are:
term1 = q_n - q**0  # q^n - 1
term2 = q_n - q**1  # q^n - q
term3 = q_n - q**2  # q^n - q^2

# Step 3: Compute the product of these terms.
# This gives the order of the general linear group GL(n, q), which is
# the number of internal adjunctions.
result = term1 * term2 * term3

# Step 4: Print the final equation and the answer.
# Remember to output each number in the final equation.
print(f"The number of internal adjunctions is the order of GL({n}, {q}).")
print(f"|GL({n}, {q})| = ({q}^{n} - {q}^0) * ({q}^{n} - {q}^1) * ({q}^{n} - {q}^2)")
print(f"|GL({n}, {q})| = ({q_n} - {q**0}) * ({q_n} - {q**1}) * ({q_n} - {q**2})")
print(f"|GL({n}, {q})| = {term1} * {term2} * {term3}")
print(f"The total number of adjunctions is: {result}")
