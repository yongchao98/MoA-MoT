import math

# Set the value for q
q = 997

# The formula for the number of involutions in PSU(4, q^2) for odd q is:
# N = (q^4 * (q^2 + 1) * (q^2 - q + 1)) / 2

# Calculate each term in the formula
term_q4 = q**4
term_q2_plus_1 = q**2 + 1
term_q2_minus_q_plus_1 = q**2 - q + 1

# Calculate the numerator
numerator = term_q4 * term_q2_plus_1 * term_q2_minus_q_plus_1

# Calculate the final result using integer division
result = numerator // 2

# Print the step-by-step calculation with the computed numbers
print(f"The number of involutions in PSU(4, {q}) is given by the formula:")
print("N = (q^4 * (q^2 + 1) * (q^2 - q + 1)) / 2")
print(f"\nSubstituting q = {q}:")
print(f"N = ({term_q4} * ({q**2} + 1) * ({q**2} - {q} + 1)) / 2")
print(f"N = ({term_q4} * {term_q2_plus_1} * {term_q2_minus_q_plus_1}) / 2")
print(f"N = {numerator} / 2")
print(f"\nNumber of involutions in PSU(4, 997): {result}")