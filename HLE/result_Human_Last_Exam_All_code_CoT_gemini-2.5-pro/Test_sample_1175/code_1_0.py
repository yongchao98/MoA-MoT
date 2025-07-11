import math

# Based on the derivation, the final expression to calculate is (3/2) * 10**(10/3) + 37/4

# Define the numbers from the final simplified equation
coeff1 = 3/2
base = 10
exponent = 10/3
const_term = 37/4

# Calculate the first term
term1 = coeff1 * (base ** exponent)

# The second term is the constant
term2 = const_term

# Calculate the final result
final_value = term1 + term2

# Output the equation and the values of its components
print("The final simplified equation is: (3/2) * 10**(10/3) + 37/4")
print(f"The first number in the equation is: {coeff1}")
print(f"The second number in the equation is: {base}")
print(f"The third number in the equation is: {exponent:.4f}")
print(f"The fourth number in the equation is: {const_term}")
print(f"The value of the final expression is: {final_value}")
