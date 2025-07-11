import math

# Step 1: Define the components to calculate p_0 = y'(0)
# The formula is p_0 = (base) ^ exponent
base = 9 / 4
exponent = 2 / 5

# Step 2: Calculate p_0 = y'(0)
p_0 = base ** exponent
print(f"The equation for the slope at x=0 is y'(0) = ({base}) ^ ({exponent})")
print(f"The value of the slope at x=0 is y'(0) = {p_0}")
print("-" * 20)

# Step 3: Define the components to calculate y(0)
# The formula is y(0) = (p_0 ^ 4) / 3
numerator = p_0 ** 4
denominator = 3

# Step 4: Calculate y(0)
y_0 = numerator / denominator
print(f"The equation for the deflection at x=0 is y(0) = (y'(0))^4 / {denominator}")
print(f"y(0) = ({p_0})^4 / {denominator}")
print(f"y(0) = {numerator} / {denominator}")
print("-" * 20)
print(f"The final deflection at x=0 is y(0) = {y_0}")
