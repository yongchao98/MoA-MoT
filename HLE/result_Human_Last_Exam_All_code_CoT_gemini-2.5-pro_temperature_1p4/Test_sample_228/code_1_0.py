import numpy as np

def create_pose_matrix(R, t):
    """Creates a 4x4 homogeneous transformation matrix (camera-to-world) from R and t."""
    T = np.eye(4)
    T[:3, :3] = R
    T[:3, 3] = t.flatten()
    return T

# Step 1: Define a 3D point in a common "world" coordinate system.
# This is the ground-truth point we are trying to reconstruct.
p_world = np.array([[5.0], [2.0], [10.0]])
p_world_h = np.vstack([p_world, [1.0]]) # Homogeneous coordinates for 4x4 matrix multiplication

# Step 2: Define the poses of two cameras relative to the world frame.
# A camera pose is the transformation that brings points from the camera's frame to the world frame.

# Camera 1: At world coordinates (-1, 0, 0)
T1_c2w = create_pose_matrix(np.eye(3), np.array([[-1.0], [0.0], [0.0]]))

# Camera 2: At world coordinates (1, 0, 0)
T2_c2w = create_pose_matrix(np.eye(3), np.array([[1.0], [0.0], [0.0]]))

# Step 3: Simulate triangulation.
# A real algorithm would use Plücker lines derived from the poses and image data.
# The crucial part is that the algorithm operates in a common frame and thus returns a point in that frame.
def triangulate_in_world_frame(pose1, pose2):
    """Dummy function to represent a triangulation algorithm."""
    print("Triangulation process running...")
    print("Geometric calculations are performed in the common 'world' frame.")
    # The result of the triangulation is a point in the world frame.
    # We use the known world point to simulate a perfect reconstruction.
    reconstructed_point_world = np.array([[5.0], [2.0], [10.0]])
    return reconstructed_point_world

# The result from triangulation is in the world frame.
p_triangulated_world = triangulate_in_world_frame(T1_c2w, T2_c2w)
p_triangulated_world_h = np.vstack([p_triangulated_world, [1.0]])


# Step 4: The result is in the world frame. To express it in a camera's frame, a transformation is needed.
# We will transform the world point into Camera 1's reference frame.
# This requires the world-to-camera transformation, which is the inverse of the camera's pose (camera-to-world).
T1_w2c = np.linalg.inv(T1_c2w)

# Apply the transformation to get the point's coordinates relative to Camera 1.
p_cam1_h = T1_w2c @ p_triangulated_world_h

print("\n--- TRANSFORMATION ---")
print("The triangulation algorithm yields a point in the world frame:")
print(f"P_world (homogeneous) = \n{p_triangulated_world_h}\n")
print("To find the point in Camera 1's frame, a transformation is required:")

print("\nEquation: P_camera1 = T_world_to_camera1 * P_world\n")
print("With the following values:")
print(f"T_world_to_camera1 = \n{np.round(T1_w2c, 2)}\n")
print("multiplied by\n")
print(f"P_world (homogeneous) = \n{p_triangulated_world_h}\n")

print("equals\n")
print(f"P_camera1 (homogeneous) = \n{np.round(p_cam1_h, 2)}\n")

print("\nThe final coordinates in Camera 1's frame are [6.0, 2.0, 10.0].")
print("This was obtained via a transformation, not directly from triangulation.")
