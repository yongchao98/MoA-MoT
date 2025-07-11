import math

# Given values
alpha = 10**16
# R = ln(100/99)
R = math.log(100/99)

# Calculate e^R and e^(2R)
exp_R = math.exp(R)
exp_2R = exp_R**2

# Calculate the numerator and denominator for t_0^2
numerator = 2 * alpha
denominator = exp_2R - 1

# Calculate t_0^2 and t_0
t0_squared = numerator / denominator
t0 = math.sqrt(t0_squared)

# Output the final equation and the result
# The problem asks to output the numbers in the final equation.
# The equation for t_0 is t_0 = sqrt(2 * alpha / (e^(2*R) - 1))
# Substituting the values: t_0 = sqrt(2 * 10^16 / ((100/99)^2 - 1))
# (100/99)^2 - 1 = 10000/9801 - 1 = 199/9801
# So, t_0 = sqrt(2 * 10^16 / (199/9801)) = sqrt(2 * 10^16 * 9801 / 199)

print("The final equation for t0 is:")
print(f"t0 = sqrt(2 * {alpha} / (e^(2*({R:.8f})) - 1))")
print(f"t0 = sqrt(2 * {alpha} / (({exp_R:.8f})^2 - 1))")
print(f"t0 = sqrt({numerator} / {denominator:.8f})")
print(f"t0 = sqrt({t0_squared})")
print("\nThe calculated positive value of t0 is:")
print(t0)