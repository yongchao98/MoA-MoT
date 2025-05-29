import numpy as np

# Define the points
C = np.array([6, 6 * np.sqrt(3), 0])
X = np.array([3, 3 * np.sqrt(3), 0])
Y = np.array([9, 3 * np.sqrt(3), 0])
Z = np.array([6, 6 * np.sqrt(3), -8])

# Function to calculate the area of a 2D triangle
def area_2d_triangle(p1, p2, p3):
    return 0.5 * abs(p1[0] * (p2[1] - p3[1]) + p2[0] * (p3[1] - p1[1]) + p3[0] * (p1[1] - p2[1]))

# Function to calculate the area of a 3D triangle
def area_3d_triangle(p1, p2, p3):
    v = p2 - p1
    w = p3 - p1
    return 0.5 * np.linalg.norm(np.cross(v, w))

# Calculate areas
area_CXY = area_2d_triangle(C, X, Y)
area_CYZ = area_3d_triangle(C, Y, Z)
area_CXZ = area_3d_triangle(C, X, Z)

# Total surface area
total_surface_area = area_CXY + area_CYZ + area_CXZ

print(total_surface_area)