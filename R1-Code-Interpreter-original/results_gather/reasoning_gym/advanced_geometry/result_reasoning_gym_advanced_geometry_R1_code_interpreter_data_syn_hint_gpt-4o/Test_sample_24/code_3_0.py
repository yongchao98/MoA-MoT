import sympy as sp

# Define the coordinates of the points
A = (-10, 5)
B = (10, 6)
C = (-1, 7)

# Define vectors BA and BC
BA = (A[0] - B[0], A[1] - B[1])
BC = (C[0] - B[0], C[1] - B[1])

# Calculate the dot product of BA and BC
dot_product = BA[0] * BC[0] + BA[1] * BC[1]

# Calculate the magnitudes of BA and BC
magnitude_BA = sp.sqrt(BA[0]**2 + BA[1]**2)
magnitude_BC = sp.sqrt(BC[0]**2 + BC[1]**2)

# Calculate the cosine of the angle at B
cos_angle_B = dot_product / (magnitude_BA * magnitude_BC)

# Calculate the angle in radians
angle_B_rad = sp.acos(cos_angle_B)

# Convert the angle to degrees
angle_B_deg = sp.deg(angle_B_rad)

# Evaluate and round the result to 3 decimal places
angle_B_deg_rounded = angle_B_deg.evalf(3)

print(angle_B_deg_rounded)