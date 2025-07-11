# The prime field characteristic is q
q = 997

# The formula for the number of involutions in PSU(4, q) for q = 1 (mod 4)
# is (q^4 * (q^2 - q + 1) * (q^2 + 1)) / 2.

# Calculate the components of the formula
q_squared = q**2
q_fourth = q**4
term1 = q_squared - q + 1
term2 = q_squared + 1

# Calculate the numerator
numerator = q_fourth * term1 * term2

# Calculate the final result using integer division
num_involutions = numerator // 2

# Print the equation with the computed values
print(f"The number of involutions is calculated by the formula:")
print(f"N = (q^4 * (q^2 - q + 1) * (q^2 + 1)) / 2")
print(f"For q = {q}:")
print(f"N = ({q_fourth} * {term1} * {term2}) / 2")
print(f"N = {numerator} / 2")
print(f"N = {num_involutions}")
