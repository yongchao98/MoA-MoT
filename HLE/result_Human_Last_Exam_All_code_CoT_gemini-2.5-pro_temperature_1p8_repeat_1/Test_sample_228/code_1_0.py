import numpy as np

def transform_point(R, t, p_world):
    """Transforms a point from world coordinates to camera coordinates."""
    # Ensure t is a column vector
    if t.ndim == 1:
        t = t.reshape(3, 1)
    # The transformation is X_cam = R * X_world + t
    # For homogenous coordinates, this is part of the view matrix
    # Here we do it manually for a 3D point
    p_cam = R @ p_world.reshape(3,1) + t
    return p_cam.flatten()

# --- 1. Define World and Camera Parameters ---
# Let's define a 3D point in the world reference frame.
P_world = np.array([2.0, 3.0, 10.0])

# Camera 1 is at the origin of the world frame.
R1 = np.eye(3)
t1 = np.zeros(3)
P1 = np.hstack([R1, t1.reshape(3, 1)]) # Projection matrix for Cam1

# Camera 2 is rotated and translated relative to Camera 1 (the world frame).
# Rotate 10 degrees around the Y axis.
theta = np.deg2rad(10)
R2 = np.array([
    [np.cos(theta), 0, np.sin(theta)],
    [0, 1, 0],
    [-np.sin(theta), 0, np.cos(theta)]
])
# Translate along the X axis.
t2 = np.array([-2.0, 0.0, 0.0])
P2 = np.hstack([R2, t2.reshape(3, 1)]) # Projection matrix for Cam2

# For simplicity, we assume an identity intrinsic matrix K (normalized coordinates).
K = np.eye(3)
P1 = K @ P1
P2 = K @ P2

# --- 2. Project the 3D point onto the 2D image planes ---
# Convert world point to homogeneous coordinates.
P_world_h = np.append(P_world, 1)

# Project onto image planes.
p1_h = P1 @ P_world_h
p2_h = P2 @ P_world_h

# Convert from homogeneous to 2D coordinates.
p1 = p1_h[:2] / p1_h[2]
p2 = p2_h[:2] / p2_h[2]

# --- 3. Triangulate the point from 2D projections ---
# The triangulation function takes the two projection matrices and the 2D points.
# This function works like Plucker coordinate-based methods by finding the best 3D point.
P_triangulated_h = cv2.triangulatePoints(P1, P2, p1, p2)

# Convert triangulated point from homogeneous to 3D coordinates.
P_triangulated_world = P_triangulated_h[:3] / P_triangulated_h[3]

# --- 4. Verify the coordinate frames ---
print(f"Original 3D Point in World Frame: {P_world}")
print(f"Triangulated 3D Point's Frame:  {P_triangulated_world.flatten()}")
print("-" * 50)
print("The triangulated point is in the World Frame (same as Camera 1's frame).")
print("To get the point in Camera 2's frame, a transformation is needed.")
print("-" * 50)

# Transform the original world point into Camera 2's frame for comparison.
P_cam2_expected = transform_point(R2, t2, P_world)
print(f"Expected Point in Camera 2 Frame: {P_cam2_expected}")

# Transform the triangulated world point into Camera 2's frame.
P_cam2_triangulated = transform_point(R2, t2, P_triangulated_world)
print(f"Triangulated Point in Camera 2 Frame (after transformation): {P_cam2_triangulated}")

print("\nConclusion: The triangulation result is in the world frame.")
print("A transformation (using rotation R2 and translation t2) was required to express it in the Camera 2 reference frame.")