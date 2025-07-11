import math

# Step 1: Define the data from the problem description and map analysis.
# Heights of the three locations in meters.
h_x = 120
h_y = 80
h_z = 140

# We measure the pixel coordinates of the points from the map image.
# Let's set a cartesian coordinate system with the origin at the bottom-left corner.
# The y-axis points North (up) and the x-axis points East (right).
# Pixel coordinates (x, y) from bottom-left corner:
X_px = (525, 405)
Y_px = (265, 305)
Z_px = (395, 590)

# The scale bar indicates that 143 pixels correspond to 200 meters.
scale_px = 143
scale_m = 200
scale_factor = scale_m / scale_px # meters per pixel

# To simplify calculations, we set Y as the origin of our map coordinates (0,0).
# We calculate the coordinates of X and Z relative to Y in meters.
x_X_rel_m = (X_px[0] - Y_px[0]) * scale_factor
y_X_rel_m = (X_px[1] - Y_px[1]) * scale_factor
x_Z_rel_m = (Z_px[0] - Y_px[0]) * scale_factor
y_Z_rel_m = (Z_px[1] - Y_px[1]) * scale_factor

# Step 2: Form vectors on the plane.
# The vectors are (dx, dy, dz), where dz is the change in height.
# Vector from Y to X
v_yx = (x_X_rel_m, y_X_rel_m, h_x - h_y)
# Vector from Y to Z
v_yz = (x_Z_rel_m, y_Z_rel_m, h_z - h_y)

# Step 3: Calculate the normal vector N = v_yx x v_yz.
# N = (A, B, C)
A = (v_yx[1] * v_yz[2]) - (v_yx[2] * v_yz[1])
B = (v_yx[2] * v_yz[0]) - (v_yx[0] * v_yz[2])
C = (v_yx[0] * v_yz[1]) - (v_yx[1] * v_yz[0])

# Step 4: Calculate the dip.
# The gradient (slope) of the plane is m = sqrt(A^2 + B^2) / |C|.
# The dip angle is arctan(m).
gradient_numerator = math.sqrt(A**2 + B**2)
gradient_denominator = abs(C)
gradient = gradient_numerator / gradient_denominator
dip_radians = math.atan(gradient)
dip_degrees = math.degrees(dip_radians)
rounded_dip = round(dip_degrees)

print("Planar Surface Dip Calculation")
print("--------------------------------")
print(f"Heights: X={h_x}m, Y={h_y}m, Z={h_z}m")
print(f"Height difference Y to X (dz_YX): {v_yx[2]} m")
print(f"Height difference Y to Z (dz_YZ): {v_yz[2]} m")
print("\nHorizontal distance from Y to X (dx, dy): ({:.1f} m, {:.1f} m)".format(v_yx[0], v_yx[1]))
print("Horizontal distance from Y to Z (dx, dy): ({:.1f} m, {:.1f} m)".format(v_yz[0], v_yz[1]))
print("\nNormal Vector N = (A, B, C)")
print(f"A = {A:.1f}")
print(f"B = {B:.1f}")
print(f"C = {C:.1f}")
print("\nThe equation for the dip angle (δ) is: δ = arctan(gradient)")
print("gradient = sqrt(A² + B²) / |C|")
print(f"gradient = sqrt({A:.1f}² + {B:.1f}²) / |{C:.1f}|")
print(f"gradient = {gradient_numerator:.1f} / {gradient_denominator:.1f} = {gradient:.4f}")
print(f"\ndip = arctan({gradient:.4f})")
print(f"dip = {dip_degrees:.2f} degrees")
print(f"\nThe dip of the planar surface, rounded to the nearest degree, is: {rounded_dip} degrees.")
print("<<<8>>>")