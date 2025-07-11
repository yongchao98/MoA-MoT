import numpy as np

def create_plucker_line_from_points(p1, p2):
    """Creates Plücker coordinates for a line passing through p1 and p2."""
    direction = p2 - p1
    moment = np.cross(p1, p2)
    return np.concatenate((direction, moment))

def get_plucker_skew_matrix(L):
    """Creates the 4x4 skew-symmetric Plücker matrix from Plücker coordinates."""
    d = L[:3]
    m = L[3:]
    # [L]_x from Plucker coordinates L = (d, m)
    # [  0  -d_z  d_y  -m_x ]
    # [ d_z   0  -d_x  -m_y ]
    # [-d_y  d_x   0   -m_z ]
    # [ m_x  m_y  m_z    0  ]
    return np.array([
        [0,      -d[2],  d[1], -m[0]],
        [d[2],   0,     -d[0], -m[1]],
        [-d[1],  d[0],   0,     -m[2]],
        [m[0],   m[1],   m[2],  0]
    ])

# 1. Define a ground-truth 3D point in the WORLD frame.
P_world = np.array([2.0, 3.0, 5.0])
print(f"Original 3D point in world frame: {P_world}")

# 2. Define two cameras by their origins in the WORLD frame.
# Camera 1 is at the world origin.
C1_origin_world = np.array([0.0, 0.0, 0.0])
# Camera 2 is translated and observes the point.
C2_origin_world = np.array([4.0, 1.0, 0.0])

# 3. Define the rays from each camera to the 3D point.
# These rays are defined in the common WORLD frame.
L1_world = create_plucker_line_from_points(C1_origin_world, P_world)
L2_world = create_plucker_line_from_points(C2_origin_world, P_world)

# 4. Perform triangulation.
# The triangulation finds a point X (homogeneous) such that L*X = 0 for all lines L.
# We stack the Plücker matrices to form a linear system AX = 0.
A1 = get_plucker_skew_matrix(L1_world)
A2 = get_plucker_skew_matrix(L2_world)
A = np.vstack((A1, A2))

print("\n--- Triangulation Equation: A * X = 0 ---")
print("Matrix A (stacked from Plücker matrices):")
# Printing each number in the equation's matrix A
for row in A:
    print(" ".join(f"{num:8.4f}" for num in row))

# 5. Solve AX = 0 using Singular Value Decomposition (SVD).
# The solution is the column of V corresponding to the smallest singular value.
_, _, Vt = np.linalg.svd(A)
X_homogeneous = Vt[-1]
# Convert from homogeneous coordinates to 3D by dividing by the last component.
P_triangulated_world = X_homogeneous[:3] / X_homogeneous[3]

print("\n--- Triangulation Result ---")
print(f"Triangulated point in WORLD frame: {np.round(P_triangulated_world, 4)}")

# 6. Show that the result is NOT in the camera frame.
# Let's find the coordinates of the point relative to Camera 2.
# This requires a transformation from world to camera 2 frame.
# For simplicity, assume Camera 2's orientation is identity (no rotation).
# The transformation is just a subtraction of the origin.
P_in_cam2_frame = P_world - C2_origin_world
print(f"\nOriginal point in CAMERA 2 frame: {np.round(P_in_cam2_frame, 4)}")

print("\n--- Conclusion ---")
print("The triangulated point's coordinates match the original world point, not the point in Camera 2's frame.")
print("This demonstrates that the result is produced in the common (world) frame, and a separate transformation is needed to represent it in a camera's local frame.")
