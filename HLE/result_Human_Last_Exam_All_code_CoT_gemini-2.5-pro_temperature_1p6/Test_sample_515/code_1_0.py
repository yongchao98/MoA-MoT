import numpy as np

# We want to calculate the probability that a 2D random walk, conditioned to avoid the origin,
# starting from (0,1), reaches the neighbors of (3600,0).
# The formula for the hitting probability of a distant point z starting from x0 is
# approximated by a(x0) / a(z), where a(x) is the potential kernel.

# The asymptotic form of the potential kernel a(x) for SRW on Z^2 is:
# a(x) ~ (2/pi) * ln|x| + C, where C is a constant.
# C = (2*gamma + ln(8)) / pi, where gamma is the Euler-Mascheroni constant.
# The target set is the four neighbors of (3600,0). We approximate this by the
# center point z = (3600,0).

# Parameters
R = 3600.0  # The distance of the target point z from the origin
x0_norm = 1.0      # The distance of the starting point x0 from the origin

# Euler-Mascheroni constant
gamma = 0.5772156649

# Calculate the constant C
C = (2 * gamma + np.log(8)) / np.pi

# The potential at the starting point x0=(0,1), |x0|=1.
# a(x0) ~ (2/pi) * ln(1) + C = C
a_x0 = C

# The potential at the target point z=(3600,0), |z|=3600.
# a(z) ~ (2/pi) * ln(3600) + C
a_z = (2 / np.pi) * np.log(R) + C

# The probability is the ratio of these potentials.
prob = a_x0 / a_z

# We need to output the equation with the calculated numbers
# P = a(x0) / a(z) = C / ( (2/pi)*ln(R) + C )
print(f"The potential kernel constant C is approximately {C:.4f}")
print(f"The potential at the starting point x0=(0,1) is a(x0) = C = {a_x0:.4f}")
print(f"The potential at the target point z=(3600,0) is a(z) = (2/pi)*ln(3600) + C = {a_z:.4f}")
print(f"The probability is the ratio a(x0) / a(z).")
print(f"P = {a_x0:.4f} / {a_z:.4f} = {prob:.4f}")

# Approximate answer with two significant digits
print(f"The approximate answer is {prob:.2g}")
# The final result is in <<<>>>
final_answer = f"{prob:.2g}"
# The final response should have each number in the equation, so we print them one by one.
# For example: print(f"2 + 2 = 4")
# We're calculating a(x0) / a(z)
print("The final calculation is:")
print(f"{a_x0} / {a_z} = {prob}")