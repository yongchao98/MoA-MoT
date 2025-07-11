import math

# Step 1 & 2: Define the geometric properties based on the plan.
# The two lines are y = x + 1 and y = -x + 5. They are perpendicular
# and intersect at P=(2, 3). The angle bisectors are x=2 and y=3.
#
# Let Circle C (radius r) have its center on one bisector (e.g., y=3).
# Let Circle D (radius 2) have its center on the other bisector (x=2).
#
# Step 3: Formulate the equations.
# The distance from a center to the vertex P is d = radius * sqrt(2).
# distance_C_P_sq = (r * sqrt(2))^2 = 2 * r^2
# distance_D_P_sq = (2 * sqrt(2))^2 = 8
#
# By Pythagorean theorem, distance(C, D)^2 = distance(C, P)^2 + distance(D, P)^2
# distance_C_D_sq = 2 * r^2 + 8
#
# The circles are tangent, so distance(C, D) = r + 2.
# Squaring this gives: distance(C, D)^2 = (r + 2)^2 = r^2 + 4r + 4
#
# Equating the two expressions for distance(C, D)^2 gives:
# 2 * r^2 + 8 = r^2 + 4r + 4
#
# Rearranging this gives the final quadratic equation: r^2 - 4r + 4 = 0.

# Step 4: Solve the quadratic equation for r and find r^2.
# The equation is in the form ar^2 + br + c = 0
a = 1
b = -4
c = 4

print(f"Based on the geometric constraints, we derive a quadratic equation for the radius r.")
# The prompt asks to output each number in the final equation.
print(f"The equation is: {a}r^2 + ({b})r + {c} = 0")

# Calculate the discriminant
discriminant = b**2 - 4*a*c

# Check if a solution exists
if discriminant < 0:
    print("The equation has no real solutions for r.")
elif discriminant == 0:
    # There is exactly one real solution
    r = -b / (2*a)
    print(f"\nSolving the equation gives a single value for the radius: r = {r}")
    
    # Calculate r squared
    r_squared = r**2
    print(f"The value of r^2 is: {r_squared}")
else:
    # There are two real solutions
    r1 = (-b - math.sqrt(discriminant)) / (2*a)
    r2 = (-b + math.sqrt(discriminant)) / (2*a)
    print(f"\nThe equation has two possible values for the radius: r = {r1} and r = {r2}")
    # The geometric context requires a single positive r, so this case wouldn't apply here,
    # but we include it for completeness of the solver.