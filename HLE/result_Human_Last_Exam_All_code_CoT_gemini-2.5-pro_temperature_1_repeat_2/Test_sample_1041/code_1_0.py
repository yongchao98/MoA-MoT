import math

# Problem parameters
R = 1000.0
z0 = (0, 300)
a1 = (0, 0)
a2 = (2, 0)

# Constants
# Euler-Mascheroni constant
gamma = 0.57721566490153286060651209008240243104215933593992

# Calculate distances
dist_z0_a1 = math.sqrt((z0[0] - a1[0])**2 + (z0[1] - a1[1])**2)
dist_z0_a2 = math.sqrt((z0[0] - a2[0])**2 + (z0[1] - a2[1])**2)
dist_a1_a2 = math.sqrt((a1[0] - a2[0])**2 + (a1[1] - a2[1])**2)

# Numerator of the probability formula
# The term is proportional to G_D(z0, a1) + G_D(z0, a2)
# We use the approximation G_D(x,y) ~ ln(R/|x-y|) for large |x-y|
num1 = math.log(R / dist_z0_a1)
num2 = math.log(R / dist_z0_a2)
numerator = num1 + num2

# Denominator of the probability formula
# The term is proportional to G_D(a1, a1) + G_D(a1, a2)
# G_D(x,x) ~ ln(R) + (pi/2)*gamma_SRW = ln(R) + gamma + ln(8)
# G_D(a1,a2) ~ ln(R) + (pi/2)*gamma_SRW - (pi/2)*a(a2-a1)
# a(2,0) = 8/(3*pi), so (pi/2)*a(2,0) = 4/3
pi_half_gamma_srw = gamma + math.log(8)
den1 = math.log(R) + pi_half_gamma_srw
den2 = math.log(R) + pi_half_gamma_srw - 4.0/3.0
denominator = den1 + den2

# Calculate the final probability
probability = numerator / denominator

# Print the equation with numerical values
# Using f-strings to format the output clearly
equation_str = (
    f"P = (log({R}/{dist_z0_a1}) + log({R}/{dist_z0_a2:.6f})) / "
    f"((log({R}) + {gamma:.6f} + log(8)) + "
    f"(log({R}) + {gamma:.6f} + log(8) - 4/3))"
)
print("The formula used for the probability is:")
print(equation_str)
print("\nWhich evaluates to:")
evaluation_str = (
    f"P = ({num1:.6f} + {num2:.6f}) / ({den1:.6f} + {den2:.6f}) "
    f"= {numerator:.6f} / {denominator:.6f}"
)
print(evaluation_str)


print(f"\nThe final probability is approximately: {probability:.3f}")
