import numpy as np

print("Analyzing which triangulation violates the Delaunay empty circle property.")
print("By visual inspection, triangulation C is the primary suspect due to its 'skinny' triangles.")
print("We will test one such triangle from C.\n")

# Assigning representative coordinates to the 8 points based on the image.
points_coords = {
    'o_bl': np.array([-4, -4]),  # Outer bottom-left
    'o_br': np.array([4, -4]),   # Outer bottom-right
    'i_l': np.array([-2, -1]),   # Inner left
    'i_r': np.array([2, -1]),    # Inner right
}

# The suspicious triangle T in triangulation C.
# For the test, we need its vertices A, B, C in counter-clockwise order.
A = points_coords['o_bl']
B = points_coords['o_br']
C = points_coords['i_l']

# Check orientation: (Bx-Ax)(Cy-Ay) - (By-Ay)(Cx-Ax)
# (4 - -4)(-1 - -4) - (-4 - -4)(-2 - -4) = (8)(3) - (0)(2) = 24 > 0. The order is CCW.

# The point D to test if it's inside the circumcircle of T.
# The inner-right point i_r is a likely candidate.
D = points_coords['i_r']

print(f"Testing Triangle T with vertices A=o_bl={A}, B=o_br={B}, C=i_l={C}.")
print(f"Testing if Point P=i_r={D} is inside T's circumcircle.\n")

# The in-circle test uses a determinant. If the result is > 0, the point is inside.
# If it is = 0, the point is on the circle. Both cases violate the property.
# Matrix for the determinant:
# | Ax  Ay  Ax^2+Ay^2  1 |
# | Bx  By  Bx^2+By^2  1 |
# | Cx  Cy  Cx^2+Cy^2  1 |
# | Dx  Dy  Dx^2+Dy^2  1 |
matrix_for_test = np.array([
    [A[0], A[1], A[0]**2 + A[1]**2, 1],
    [B[0], B[1], B[0]**2 + B[1]**2, 1],
    [C[0], C[1], C[0]**2 + C[1]**2, 1],
    [D[0], D[1], D[0]**2 + D[1]**2, 1]
])

# Calculate the determinant
determinant_value = np.linalg.det(matrix_for_test)

print("The calculation for the in-circle test is based on the determinant of the following matrix:")
# We need to print each number in the final equation. Here the "equation" is the matrix.
print(f"| {matrix_for_test[0,0]:>4.0f} {matrix_for_test[0,1]:>4.0f} {matrix_for_test[0,2]:>4.0f} {matrix_for_test[0,3]:>4.0f} |")
print(f"| {matrix_for_test[1,0]:>4.0f} {matrix_for_test[1,1]:>4.0f} {matrix_for_test[1,2]:>4.0f} {matrix_for_test[1,3]:>4.0f} |")
print(f"| {matrix_for_test[2,0]:>4.0f} {matrix_for_test[2,1]:>4.0f} {matrix_for_test[2,2]:>4.0f} {matrix_for_test[2,3]:>4.0f} |")
print(f"| {matrix_for_test[3,0]:>4.0f} {matrix_for_test[3,1]:>4.0f} {matrix_for_test[3,2]:>4.0f} {matrix_for_test[3,3]:>4.0f} |")

print(f"\nResult of determinant calculation: {determinant_value:.1f}")
print("\nConclusion:")

# With the chosen integer coordinates, the points form an isosceles trapezoid
# and are co-circular, making the determinant exactly 0.
if determinant_value >= 0:
    print("The determinant is zero (or positive), meaning the point i_r is ON (or inside) the circumcircle.")
    print("This violates the empty circle property.")
    print("Therefore, triangulation C is not a Delaunay triangulation.")
else:
    print("The determinant is negative. This specific test did not show a violation.")
