import numpy as np
from scipy.spatial.transform import Rotation as R

def cross_product_matrix(v):
    """Returns the skew-symmetric cross-product matrix for a 3D vector v."""
    return np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

# Step 1: Define a 3D point in the world reference frame.
X_world = np.array([2.0, 3.0, 5.0])

print("--- Setup ---")
print(f"Original 3D point in WORLD frame (X_world): {X_world}\n")

# Step 2: Define two cameras. Camera 1 is at the world origin.
# Camera 2 is translated and rotated.
R1 = np.identity(3)
C1 = np.array([0.0, 0.0, 0.0])
t1 = -R1 @ C1

r2_rot = R.from_euler('y', 30, degrees=True)
R2 = r2_rot.as_matrix()
C2 = np.array([10.0, 0.0, 0.0])
t2 = -R2 @ C2

print("--- Camera Poses (in World Frame) ---")
print(f"Camera 1 Center (C1): {C1}")
print(f"Camera 2 Center (C2): {C2}\n")


# Step 3: Construct the two projection rays (lines) in the WORLD frame.
# A Plücker line is defined by a direction vector (d) and a moment vector (m = p x d),
# where p is any point on the line (we use the camera centers).
print("--- Ray Construction in World Frame ---")
# Ray 1
d1_world = (X_world - C1) / np.linalg.norm(X_world - C1)
m1_world = np.cross(C1, d1_world)
print(f"Ray 1 (L1) Plücker Coords [d1, m1] in world frame: {np.round(np.concatenate((d1_world, m1_world)), 3)}")

# Ray 2
d2_world = (X_world - C2) / np.linalg.norm(X_world - C2)
m2_world = np.cross(C2, d2_world)
print(f"Ray 2 (L2) Plücker Coords [d2, m2] in world frame: {np.round(np.concatenate((d2_world, m2_world)), 3)}\n")


# Step 4: Triangulate the point by finding the intersection of the two world-frame lines.
# The result will be in the same frame as the lines: the world frame.
print("--- Triangulation ---")
# We solve the system [d1_x]X = m1 and [d2_x]X = m2 using least squares.
A = np.vstack([cross_product_matrix(d1_world), cross_product_matrix(d2_world)])
b = np.concatenate([m1_world, m2_world])
X_reconstructed_world = np.linalg.pinv(A) @ b

print(f"Reconstructed 3D point: {np.round(X_reconstructed_world, 3)}")
print("This point is in the WORLD reference frame.\n")


# Step 5: To get the point in a camera's frame, a transformation is needed.
print("--- Transformation to Camera 2 Frame ---")
print("To find the point's coordinates relative to Camera 2, we apply its pose transformation:")
print("Equation: X_camera2 = R2 * X_world + t2")

# Apply the transformation to get the point in Camera 2's reference frame
X_in_cam2_frame = R2 @ X_reconstructed_world + t2

# Print the final equation with all the numbers
print("\nFinal Calculation:")
print("X_camera2 = ")
print(f"{np.round(R2, 4)}\n    (R2)\n")
print("          *\n")
print(f"{np.round(X_reconstructed_world, 4)}\n    (X_world)\n")
print("          +\n")
print(f"{np.round(t2, 4)}\n    (t2)\n")
print("---------------------------------")
print(f"Result in Camera 2 Frame = {np.round(X_in_cam2_frame, 4)}")
print("\nConclusion: The transformation was a necessary, separate step.")
