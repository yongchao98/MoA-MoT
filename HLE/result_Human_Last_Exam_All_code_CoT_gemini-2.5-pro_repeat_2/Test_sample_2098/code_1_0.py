import sympy as sp

# Define the variable and the function for the robotic arm's path
x = sp.Symbol('x')
# Based on analysis, we use a plausible intended form of y1(x)
y1 = 16/x**2 - 32/x**3

# Define the integration limits
lower_limit = 2
upper_limit = 4

# Calculate the integral
integral_val = sp.integrate(y1, (x, lower_limit, upper_limit))

# The integral is I = [-16/x + 16/x^2] from 2 to 4
# I = (-16/4 + 16/16) - (-16/2 + 16/4)
# I = (-4 + 1) - (-8 + 4)
# I = -3 - (-4) = 1

# We print the components of the final calculation
term1_at_4 = -16/4
term2_at_4 = 16/(4**2)
val_at_4 = term1_at_4 + term2_at_4

term1_at_2 = -16/2
term2_at_2 = 16/(2**2)
val_at_2 = term1_at_2 + term2_at_2

result = val_at_4 - val_at_2

print(f"The integral is calculated as the definite integral of (16/x^2 - 32/x^3) from {lower_limit} to {upper_limit}.")
print(f"The indefinite integral is F(x) = -16/x + 16/x^2.")
print(f"Evaluating at the upper limit: F({upper_limit}) = -16/{upper_limit} + 16/{upper_limit}^2 = {val_at_4}")
print(f"Evaluating at the lower limit: F({lower_limit}) = -16/{lower_limit} + 16/{lower_limit}^2 = {val_at_2}")
print(f"The result of the integral is F(4) - F(2) = ({val_at_4}) - ({val_at_2}) = {result}")
