import math

# Based on the asymptotic analysis, the solution takes the form y(x) = A * x^p.
# The analysis revealed the values for p and A^3.
p = -1.00
A_cubed = -6.00

# Calculate the coefficient A
A = math.copysign(1, A_cubed) * (abs(A_cubed)**(1/3))

# Round A to two decimal places for the final expression
A_rounded = round(A, 2)

# Output the components of the derived analytical expression
print(f"The analysis suggests an approximate solution of the form y(x) = A * x^p.")
print(f"The exponent is p = {p:.2f}")
print(f"The coefficient A is approximately A = {A_rounded}")

# Display the final analytical expression as requested
numerator = A_rounded
print(f"\nThe resulting analytical expression for y(x) in the large x regime is:")
print(f"y(x) = {numerator} / x")