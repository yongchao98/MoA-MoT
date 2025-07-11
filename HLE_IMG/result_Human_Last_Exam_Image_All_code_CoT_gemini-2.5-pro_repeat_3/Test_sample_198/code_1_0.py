import math
import numpy as np

# Step 1: Define the data from the problem statement and map analysis.
# Heights of the three locations in meters.
h_X = 120.0
h_Y = 80.0
h_Z = 140.0

# Pixel coordinates measured from the image (origin at top-left).
# We assume the y-axis points up (North) and x-axis points right (East).
px_X = (543, 461)
px_Y = (254, 577)
px_Z = (452, 252)

# Pixel coordinates of the scale bar.
scale_start_px = (757, 604)
scale_end_px = (932, 604)
scale_dist_m = 200.0

# Calculate the scale factor (meters per pixel).
scale_dist_px = scale_end_px[0] - scale_start_px[0]
scale_factor = scale_dist_m / scale_dist_px

print(f"Step 1: Coordinate and Scale Determination")
print(f"Height of X: {h_X} m, Y: {h_Y} m, Z: {h_Z} m")
print(f"Scale bar is {scale_dist_px} pixels long, representing {scale_dist_m} meters.")
print(f"Scale factor: {scale_factor:.4f} m/pixel\n")


# Step 2: Establish a map coordinate system with Y at the origin (0,0) and calculate coordinates in meters.
# We set the y-axis to point North (up), so we subtract y-pixel coordinates from a larger value (like the origin's).
x_Y_m, y_Y_m = 0.0, 0.0
x_X_m = (px_X[0] - px_Y[0]) * scale_factor
y_X_m = (px_Y[1] - px_X[1]) * scale_factor # Invert y-axis for North up
x_Z_m = (px_Z[0] - px_Y[0]) * scale_factor
y_Z_m = (px_Y[1] - px_Z[1]) * scale_factor # Invert y-axis for North up

# Create 3D points.
p_X = np.array([x_X_m, y_X_m, h_X])
p_Y = np.array([x_Y_m, y_Y_m, h_Y])
p_Z = np.array([x_Z_m, y_Z_m, h_Z])

print(f"Step 2: 3D Coordinates (Y at origin)")
print(f"Point X: ({p_X[0]:.1f} m, {p_X[1]:.1f} m, {p_X[2]:.1f} m)")
print(f"Point Y: ({p_Y[0]:.1f} m, {p_Y[1]:.1f} m, {p_Y[2]:.1f} m)")
print(f"Point Z: ({p_Z[0]:.1f} m, {p_Z[1]:.1f} m, {p_Z[2]:.1f} m)\n")

# Step 3: Create two vectors on the plane and find the normal vector.
v_YX = p_X - p_Y
v_YZ = p_Z - p_Y

# Calculate the cross product to get the normal vector n = (nx, ny, nz).
n = np.cross(v_YX, v_YZ)
nx, ny, nz = n[0], n[1], n[2]

print(f"Step 3: Normal Vector Calculation")
print(f"Vector YX: <{v_YX[0]:.1f}, {v_YX[1]:.1f}, {v_YX[2]:.1f}>")
print(f"Vector YZ: <{v_YZ[0]:.1f}, {v_YZ[1]:.1f}, {v_YZ[2]:.1f}>")
print(f"Normal vector n = YX x YZ = <{nx:.1f}, {ny:.1f}, {nz:.1f}>\n")

# Step 4: Calculate the dip angle.
# tan(dip) = sqrt(nx^2 + ny^2) / |nz|
tan_dip = math.sqrt(nx**2 + ny**2) / abs(nz)
dip_rad = math.atan(tan_dip)
dip_deg = math.degrees(dip_rad)

# Round the final answer to the nearest degree.
rounded_dip = round(dip_deg)

print(f"Step 4: Dip Calculation")
print("The dip angle (δ) is calculated using the formula: δ = arctan(sqrt(nx² + ny²) / |nz|)")
print(f"Plugging in the values:")
print(f"δ = arctan(sqrt(({nx:.1f})² + ({ny:.1f})²) / |{nz:.1f}|)")
print(f"δ = arctan(sqrt({nx**2:.1f} + {ny**2:.1f}) / {abs(nz):.1f})")
print(f"δ = arctan({tan_dip:.4f})")
print(f"δ = {dip_deg:.2f} degrees\n")

# Step 5: Final Answer
print(f"Step 5: Final Result")
print(f"The calculated dip of the planar surface is {dip_deg:.2f} degrees.")
print(f"Rounding to the nearest degree, the dip is {rounded_dip} degrees.")
print(f"<<<{rounded_dip}>>>")