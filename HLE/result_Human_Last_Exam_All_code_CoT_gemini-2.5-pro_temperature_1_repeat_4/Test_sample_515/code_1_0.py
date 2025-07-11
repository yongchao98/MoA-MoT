import math

# Constants
# gamma is the Euler-Mascheroni constant
gamma = 0.5772156649
# z is the target point (3600, 0)
z_norm = 3600
# x0 is the starting point (0, 1)
x0_norm = 1

# The potential kernel a(x) is approximated by (2/pi)*ln(||x||) + C for large ||x||.
# C = (2*gamma + ln(8))/pi
# The probability is given by a(x0) / a(z).
# We approximate a(x0) with the formula at ||x0||=1, which gives a(x0) approx C, as ln(1)=0.

# Numerator of the probability is a(x0)
# We use the asymptotic formula for consistency.
# a(x0) approx (2/pi)*ln(1) + (2*gamma + ln(8))/pi = (2*gamma + ln(8))/pi
numerator_term_without_pi = 2 * gamma + math.log(8)

# Denominator of the probability is a(z)
# a(z) approx (2/pi)*ln(3600) + (2*gamma + ln(8))/pi
denominator_term_without_pi = 2 * math.log(z_norm) + (2 * gamma + math.log(8))

# The (1/pi) factor cancels out in the ratio.
probability = numerator_term_without_pi / denominator_term_without_pi

# Print the formula with the calculated values
print(f"The probability is calculated as the ratio of potential kernels a(x0)/a(z), where x0=(0,1) and z=(3600,0).")
print(f"Using the asymptotic formula a(x) ~ (2/pi)*ln(||x||) + C, the probability is C / ((2/pi)*ln(3600) + C).")
print(f"This simplifies to (2*gamma + ln(8)) / (2*ln(3600) + 2*gamma + ln(8)).")
print(f"P = ({2 * gamma:.4f} + {math.log(8):.4f}) / (2*{math.log(z_norm):.4f} + {2 * gamma:.4f} + {math.log(8):.4f})")
print(f"P = {numerator_term_without_pi:.4f} / ({2 * math.log(z_norm):.4f} + {numerator_term_without_pi:.4f})")
print(f"P = {numerator_term_without_pi:.4f} / {denominator_term_without_pi:.4f}")
print(f"P = {probability:.4f}")

# Print the final answer with two significant digits.
print(f"\nThe approximate answer with two significant digits is: {probability:.2g}")