import numpy as np

def triangulate_plucker(P1, P2, x1, x2):
    """
    Triangulates a 3D point from two 2D correspondences using Plücker coordinates.
    The calculation is performed in the world frame.
    
    Args:
        P1, P2: 3x4 projection matrices for camera 1 and 2.
        x1, x2: 2D points in homogeneous coordinates (3x1 vector).
        
    Returns:
        The triangulated 3D point in the world frame as a 3x1 vector.
    """

    # Helper function to back-project a 2D point to a 3D Plücker line in the world frame.
    def back_project_to_plucker(P, x):
        # The camera center C is the null space of the projection matrix P.
        U, S, Vt = np.linalg.svd(P)
        C = Vt[-1, :].T
        C = C / C[3]  # Dehomogenize to get 3D coordinates in world frame

        # Find another point on the ray by back-projecting the image point x.
        X_on_ray = np.linalg.pinv(P) @ x
        X_on_ray = X_on_ray / X_on_ray[3]  # Dehomogenize

        # A Plücker line L is defined by a direction vector 'd' and a moment vector 'm'.
        # Both are calculated here in the world frame.
        d = X_on_ray[:3] - C[:3]
        d = d / np.linalg.norm(d)
        m = np.cross(C[:3], d)
        return d, m

    # 1. Get the Plücker coordinates for the two rays in the world frame.
    d1, m1 = back_project_to_plucker(P1, x1)
    d2, m2 = back_project_to_plucker(P2, x2)

    # 2. Find the point of closest approach between the two skew lines.
    # This point is the triangulated 3D point in the world frame.
    # In the ideal case (no noise), the lines intersect and this gives the intersection point.
    d1_cross_d2 = np.cross(d1, d2)
    d1_cross_d2_sq_norm = np.dot(d1_cross_d2, d1_cross_d2)

    # Handle the case of parallel lines
    if d1_cross_d2_sq_norm < 1e-9:
        print("Error: Lines are parallel, triangulation is ambiguous.")
        return None

    # Calculate the closest point on each line to the other line.
    p1 = (np.dot(m1, d2) * d1 + np.dot(d1, d2) * m1 - np.dot(d1, m1) * d2) / d1_cross_d2_sq_norm
    p2 = (-np.dot(m2, d1) * d2 - np.dot(d1, d2) * m2 + np.dot(d2, m2) * d1) / d1_cross_d2_sq_norm

    # The optimal 3D point is the midpoint of the segment connecting p1 and p2.
    X_world = (p1 + p2) / 2.0
    return X_world

# --- Main Demonstration ---

# Define a true 3D point in the WORLD coordinate frame.
X_world_true = np.array([2, 3, 10, 1])  # Homogeneous coordinates

# Define two cameras.
# Camera 1 is at the origin of the world frame. Its frame is the world frame.
R1 = np.eye(3)
t1 = np.zeros((3, 1))

# Camera 2 is translated 5 units along the x-axis relative to the world frame.
R2 = np.eye(3)
t2 = np.array([[-5], [0], [0]])

# Define a simple intrinsic matrix K (same for both cameras).
K = np.array([[1000, 0, 500], [0, 1000, 500], [0, 0, 1]])

# Create the full projection matrices P = K[R|t].
P1 = K @ np.hstack((R1, t1))
P2 = K @ np.hstack((R2, t2))

# Project the 3D world point onto each camera's image plane.
x1_proj = P1 @ X_world_true
x1 = x1_proj / x1_proj[2]  # Normalize to get homogeneous 2D coordinates

x2_proj = P2 @ X_world_true
x2 = x2_proj / x2_proj[2]  # Normalize

# Perform triangulation. This function returns the point in the WORLD frame.
X_triangulated_world = triangulate_plucker(P1, P2, x1, x2)

print(f"Original 3D point in world frame: {X_world_true[:3]}")
print(f"Triangulated 3D point in world frame: {np.round(X_triangulated_world, 5)}")
print("-" * 50)
print("The result above is in the world frame, NOT a camera's local frame.")
print("To get the point in a camera's frame, a transformation is needed.")
print("-" * 50)

# To get the point's coordinates in camera 2's frame, we apply its extrinsic transformation.
# X_cam2 = R2 * X_world + t2
X_in_cam2_frame = R2 @ X_triangulated_world + t2.flatten()

# The point's coordinates in camera 1's frame are the same as the world frame in this specific setup.
X_in_cam1_frame = R1 @ X_triangulated_world + t1.flatten()

print(f"Point coordinates in Camera 1's reference frame: {np.round(X_in_cam1_frame, 5)}")
print(f"Point coordinates in Camera 2's reference frame: {np.round(X_in_cam2_frame, 5)}")
