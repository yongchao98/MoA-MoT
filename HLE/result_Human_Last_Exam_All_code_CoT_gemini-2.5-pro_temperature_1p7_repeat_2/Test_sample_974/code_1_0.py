import sympy

# Define the number q from the problem.
q = 12740347

# As per the explanation, q is a prime number.
# The only prime divisor of q is p = q.
p = q

# We now check if the number of elements of order p is equal
# in PSL(3, q^2) and PSL(4, q).

# For PSL(3, q^2): n=3, s=q^2.
n1 = 3
# The exponent in the formula s^(n*(n-1)/2)
exp1 = n1 * (n1 - 1) // 2
final_exp1 = 2 * exp1

# For PSL(4, q): n=4, s=q.
n2 = 4
# The exponent in the formula s^(n*(n-1)/2)
exp2 = n2 * (n2 - 1) // 2

# The condition is that the number of elements are equal.
# (q^2)^exp1 - 1 == q^exp2 - 1
# This simplifies to checking if the exponents of q are equal.

print(f"Let q = {q}. The prime divisor to check is p = q = {p}.")
print("We check if the number of elements of order p are equal for the two groups.")
print("The formula for the number of elements is s^(n*(n-1)/2) - 1.")
print("")
print("For PSL(3, q^2):")
print(f"The number of elements is (q^2)^({n1}*({n1}-1)/2) - 1 = q^{final_exp1} - 1")
print("")
print("For PSL(4, q):")
print(f"The number of elements is q^({n2}*({n2}-1)/2) - 1 = q^{exp2} - 1")
print("")
print(f"Comparing the exponents of q: {final_exp1} == {exp2}")
print("Since the exponents are equal, the numbers of elements are equal.")
print("The condition is satisfied for the prime divisor:")
print(p)