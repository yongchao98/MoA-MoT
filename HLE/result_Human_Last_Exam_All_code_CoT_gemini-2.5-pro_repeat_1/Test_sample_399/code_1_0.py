# The value of t for which we want to find x(t)
t = 1

# The equation for x(t) is x(t) = 8 / (-4*t^2 + 8*t - 1)
# Let's calculate the numerator and denominator for t=1

# Numerator is constant
numerator = 8

# Terms in the denominator
term1_coeff = -4
term2_coeff = 8
term3_coeff = -1

# Calculate the denominator value
denominator_val = term1_coeff * (t**2) + term2_coeff * t + term3_coeff

# Calculate the final value of x(1)
x_1 = numerator / denominator_val

# Print the final equation with the numbers plugged in
print(f"The solution to the differential equation is x(t) = 8 / (-4*t^2 + 8*t - 1).")
print(f"To find x(1), we substitute t=1 into the equation:")
print(f"x(1) = {numerator} / ({term1_coeff}*({t})^2 + {term2_coeff}*({t}) + {term3_coeff})")
print(f"x(1) = {numerator} / ({term1_coeff * (t**2)} + {term2_coeff * t} + {term3_coeff})")
print(f"x(1) = {numerator} / {denominator_val}")
print(f"The final value is: {x_1}")