# This script prints the formula for the correction factor P(n).

# Define the coefficients for the formula of P(n) of the form:
# P(n) = c0 + c1 / n^2 + c2 / n^4
c0 = 1
c1_num = 1
c1_den = 720
c2_num = -1433
c2_den = 7257600

# Print the final formula in a structured way.
print("The formula for P(n) is:")
print(f"P(n) = {c0} + ({c1_num})/({c1_den} * n**2) + ({c2_num})/({c2_den} * n**4)")

# Print the formula in a more commonly seen format.
print("\nThis can be written as:")
print(f"P(n) = 1 + {c1_num}/({c1_den}*n**2) - {abs(c2_num)}/({c2_den}*n**4)")

# Outputting each number in the final equation as per the instructions.
print("\nThe numbers in the final equation are:")
print(f"Term 1 (constant): {c0}")
print(f"Term 2 (coefficient of 1/n^2): {c1_num}/{c1_den}")
print(f"Term 3 (coefficient of 1/n^4): {c2_num}/{c2_den}")
