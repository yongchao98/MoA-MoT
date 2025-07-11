import math

# Given parameters
alpha = 10**16
# R = ln(100/99)
R = math.log(100/99)

# The solvability condition leads to the equation:
# t_0^2 * (e^R - 1) / 5 = alpha
# We solve for t_0.

# Calculate e^R - 1
eR_minus_1 = math.exp(R) - 1

# The equation is t_0^2 * (eR_minus_1 / 5) = alpha
# t_0^2 = 5 * alpha / eR_minus_1
t0_squared = 5 * alpha / eR_minus_1

# Calculate the positive value of t_0
t0 = math.sqrt(t0_squared)

# The problem asks to output each number in the final equation.
# The equation is t_0^2 * (1/99) * (1/5) = 10^16
val_eR_minus_1_num = 1
val_eR_minus_1_den = 99
val_integral = 5
val_t0_squared_factor = val_integral * val_eR_minus_1_den

print("The final equation for t_0^2 is derived from the solvability condition:")
print(f"t_0^2 * (1/{val_eR_minus_1_den}) * (1/{val_integral}) = {alpha:.0e}")
print(f"t_0^2 = {val_integral} * {val_eR_minus_1_den} * {alpha:.0e}")
print(f"t_0^2 = {val_t0_squared_factor} * {alpha:.0e}")
print(f"t_0 = sqrt({val_t0_squared_factor}) * 10^8")
print(f"The positive value of t_0 is: {t0}")
