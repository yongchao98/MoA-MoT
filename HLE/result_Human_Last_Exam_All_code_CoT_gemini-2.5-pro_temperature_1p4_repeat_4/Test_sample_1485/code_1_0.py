import numpy as np

# This script simulates the rotation of a 3D torus and renders it as ASCII art.

# --- 1. Define Parameters ---

# Torus geometry, estimated from the initial image
R = 7.0  # Major radius (center of torus to center of tube)
r = 4.0  # Minor radius (radius of the tube)

# Rotation angles in degrees for X, Y, and Z axes
rot_x, rot_y, rot_z = 140, 75, 35

# Screen (canvas) dimensions for the output
width, height = 50, 22

# --- 2. Generate Torus Points ---

# Create a grid of angles for the parametric equations
num_u = 180  # Density of points around the major circle
num_v = 90   # Density of points around the minor circle
u = np.linspace(0, 2 * np.pi, num_u)
v = np.linspace(0, 2 * np.pi, num_v)
u, v = np.meshgrid(u, v)

# Parametric equations to define the surface of the torus
x = (R + r * np.cos(v)) * np.cos(u)
y = (R + r * np.cos(v)) * np.sin(u)
z = r * np.sin(v)

# Combine all points into a single (3, N) matrix for efficient processing
points = np.stack((x.flatten(), y.flatten(), z.flatten()))

# --- 3. Create and Apply Rotation Matrix ---

# Convert degrees to radians. Use negative angles because the problem
# defines positive rotation as clockwise, while standard matrices are counter-clockwise.
ax_rad = np.deg2rad(-rot_x)
ay_rad = np.deg2rad(-rot_y)
az_rad = np.deg2rad(-rot_z)

# Standard 3D rotation matrices
Rx = np.array([
    [1, 0, 0],
    [0, np.cos(ax_rad), -np.sin(ax_rad)],
    [0, np.sin(ax_rad), np.cos(ax_rad)]
])
Ry = np.array([
    [np.cos(ay_rad), 0, np.sin(ay_rad)],
    [0, 1, 0],
    [-np.sin(ay_rad), 0, np.cos(ay_rad)]
])
Rz = np.array([
    [np.cos(az_rad), -np.sin(az_rad), 0],
    [np.sin(az_rad), np.cos(az_rad), 0],
    [0, 0, 1]
])

# Combine matrices. Rotation order is X, then Y, then Z (P' = Rz * Ry * Rx * P)
rotation_matrix = Rz @ Ry @ Rx

# Apply the combined rotation to all points of the torus
rotated_points = rotation_matrix @ points

# --- 4. Project onto 2D Screen and Z-Buffer ---

xp, yp, zp = rotated_points[0, :], rotated_points[1, :], rotated_points[2, :]

# Scale and shift points to fit the canvas dimensions
scale = width / (2.2 * (R + r))
x_img = (xp * scale + width / 2).astype(int)
y_img = (-yp * scale + height / 2).astype(int)  # Y is inverted for screen coordinates

# Initialize canvas and z-buffer (for hidden surface removal)
canvas = np.full((height, width), ' ')
zbuffer = np.full((height, width), np.inf)

# Z-buffering: only keep the point closest to the observer for each pixel
for i in range(len(x_img)):
    xi, yi = x_img[i], y_img[i]
    if 0 <= xi < width and 0 <= yi < height:
        if zp[i] < zbuffer[yi, xi]:
            zbuffer[yi, xi] = zp[i]

# --- 5. Render Final ASCII Art ---

# Define shades: from nearest ('░') to farthest ('█')
shades = ['░', '▒', '▓', '█']

# Find the depth range of visible points to normalize for shading
visible_z = zbuffer[zbuffer != np.inf]
if len(visible_z) > 0:
    min_z, max_z = np.min(visible_z), np.max(visible_z)
    z_range = max_z - min_z if (max_z - min_z) > 0 else 1
else:
    min_z, z_range = 0, 1

# Fill the canvas with shaded characters based on depth
for r_idx in range(height):
    for c_idx in range(width):
        if zbuffer[r_idx, c_idx] != np.inf:
            norm_z = (zbuffer[r_idx, c_idx] - min_z) / z_range
            shade_index = min(int(norm_z * len(shades)), len(shades) - 1)
            canvas[r_idx, c_idx] = shades[shade_index]

# --- 6. Print the Output ---

print(f"Final view for rotation (X, Y, Z) = ({rot_x}, {rot_y}, {rot_z}):")
for row in canvas:
    print("".join(row))
