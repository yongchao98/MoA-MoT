import math

# The problem is to evaluate the sum S = sum_{n=2 to infinity} (-1)^n * S_n / n
# As derived, the value is S = (1/2)*ln(2)^2 + 1 - pi^2/12

# Calculate the components of the final equation
term1 = 0.5 * (math.log(2))**2
term2 = 1.0
term3 = -(math.pi**2) / 12

# The final equation is S = term1 + term2 + term3
final_value = term1 + term2 + term3

# Print each number in the final equation
print("The value of the first term (1/2 * ln(2)^2) is:")
print(term1)

print("\nThe value of the second term (1) is:")
print(term2)

print("\nThe value of the third term (-pi^2 / 12) is:")
print(term3)

print("\nThe final value of the sum is:")
print(final_value)
