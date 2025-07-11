import numpy as np

def calculate_angle(p1, p2, p3):
    """Calculates the angle p1-p2-p3 at vertex p2 in degrees."""
    v1 = p1 - p2
    v2 = p3 - p2
    dot_product = np.dot(v1, v2)
    mag_v1 = np.linalg.norm(v1)
    mag_v2 = np.linalg.norm(v2)
    cosine_angle = dot_product / (mag_v1 * mag_v2)
    # Clip to handle potential floating point inaccuracies leading to values > 1 or < -1
    cosine_angle = np.clip(cosine_angle, -1.0, 1.0)
    angle_rad = np.arccos(cosine_angle)
    return np.degrees(angle_rad)

# Define coordinates for the points based on a regular pentagon of radius 100,
# and two inner points.
# v1: top, v2: top-right, v3: bottom-right
# i2: inner-right point
R = 100
v1 = np.array([R * np.cos(np.deg2rad(90)), R * np.sin(np.deg2rad(90))])
v2 = np.array([R * np.cos(np.deg2rad(18)), R * np.sin(np.deg2rad(18))])
v3 = np.array([R * np.cos(np.deg2rad(-54)), R * np.sin(np.deg2rad(-54))])
i2 = np.array([30, 0])

# In triangulation B, consider the convex quadrilateral v1-v2-v3-i2.
# It is divided by the diagonal edge v1-v3.
# The two triangles are (v1, v2, v3) and (v1, i2, v3).
# We check the sum of the angles opposite the shared edge v1-v3.

# Angle at vertex v2 in the quadrilateral
angle1 = calculate_angle(v1, v2, v3)

# Angle at vertex i2 in the quadrilateral
angle2 = calculate_angle(v1, i2, v3)

total_angle = angle1 + angle2

print("Analysis of Triangulation B:")
print("Considering the quadrilateral formed by vertices v1, v2, v3, and i2.")
print("The triangulation uses the diagonal edge v1-v3.")
print("To check if this edge is Delaunay, we sum the opposite angles.")
print(f"The angle at vertex v2 (angle v1-v2-v3) is: {angle1:.2f} degrees.")
print(f"The angle at vertex i2 (angle v1-i2-v3) is: {angle2:.2f} degrees.")
print(f"The sum of these opposite angles is: {angle1:.2f} + {angle2:.2f} = {total_angle:.2f} degrees.")
print(f"Since the sum ({total_angle:.2f}) is greater than 180 degrees, the edge v1-v3 is an 'illegal' edge.")
print("Therefore, triangulation B violates the empty circle property of Delaunay triangulations.")
