import math

# Based on the derivation, the area of the triangle T(t) is a function of time.
# The formula is determined by the movement of the triangle's vertices along
# the sides of the hexagon. The rotation of the hexagon does not affect the area.

# The final formula for the area A(t) involves several constants derived from
# the geometry and motion of the vertices. The Python script below constructs
# and prints this formula, ensuring every number in the equation is explicitly shown.

# Define the numerical constants from the derived formula:
# A(t) = (c1_num/c1_den)*sqrt(3) + (c2_num/c2_den)*sqrt(3) * (term1 - abs(term2 - (t % term3_mod)))**2
c1_num = 225
c1_den = 4
c2_num = 3
c2_den = 4
term1 = 5
term2 = 5
term3_mod = 10

# Print the final function for the area A(t).
# The f-string is used to embed the numerical constants directly into the output string.
print("The area of the triangle T(t) as a function of time t is given by the formula:")
print(
    f"A(t) = ({c1_num}/{c1_den}) * sqrt(3) + ({c2_num}/{c2_den}) * sqrt(3) * "
    f"({term1} - abs({term2} - (t % {term3_mod})))**2"
)

# To demonstrate the periodic nature of the area, we can define a function
# to calculate the area at specific times.
def calculate_area(t):
    """Calculates the area of triangle T(t) for a given time t."""
    delta_t_abs = term1 - abs(term2 - (t % term3_mod))
    area = (c1_num / c1_den) * math.sqrt(3) + (c2_num / c2_den) * math.sqrt(3) * (delta_t_abs**2)
    return area

print("\nFor example, we can calculate the area at different times:")
print(f"Area at t = 0.0 seconds: {calculate_area(0.0):.4f}")
print(f"Area at t = 5.0 seconds (max area): {calculate_area(5.0):.4f}")
print(f"Area at t = 10.0 seconds (min area): {calculate_area(10.0):.4f}")
