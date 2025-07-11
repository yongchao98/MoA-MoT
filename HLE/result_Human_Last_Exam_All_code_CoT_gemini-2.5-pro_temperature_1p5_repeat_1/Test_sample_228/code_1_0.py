import numpy as np

def skew_symmetric(v):
    """Creates a skew-symmetric matrix from a 3-element vector."""
    return np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

def triangulate_point(ray1_origin, ray1_dir, ray2_origin, ray2_dir):
    """
    Finds the point of closest approach between two rays in 3D space.
    This point is the triangulated 3D point.
    """
    # System of equations to solve for the closest point
    # Based on the principle that the line connecting the closest points
    # is perpendicular to both ray directions.
    ray1_dir = ray1_dir / np.linalg.norm(ray1_dir)
    ray2_dir = ray2_dir / np.linalg.norm(ray2_dir)

    A = np.array([
        [np.dot(ray1_dir, ray1_dir), -np.dot(ray1_dir, ray2_dir)],
        [np.dot(ray2_dir, ray1_dir), -np.dot(ray2_dir, ray2_dir)]
    ])
    
    b = np.array([
        np.dot(ray1_dir, ray2_origin - ray1_origin),
        np.dot(ray2_dir, ray2_origin - ray1_origin)
    ])
    
    try:
        # t1, t2 are the distances along each ray to the points of closest approach
        t1, t2 = np.linalg.solve(A, b)
        
        # The two closest points on each ray
        p1_closest = ray1_origin + t1 * ray1_dir
        p2_closest = ray2_origin + t2 * ray2_dir
        
        # The triangulated point is the midpoint between these two points
        triangulated_point = (p1_closest + p2_closest) / 2
        return triangulated_point
        
    except np.linalg.LinAlgError:
        print("Error: The rays are parallel, cannot solve for a unique intersection.")
        return None

# 1. Define a 3D point in the world coordinate system
P_world = np.array([2.0, 3.0, 10.0])

# 2. Define two cameras
# Camera 1 is at the origin of the world frame
R1_world = np.identity(3)
T1_world = np.array([0., 0., 0.])
C1_world = -R1_world.T @ T1_world # Camera center in world frame

# Camera 2 is translated and rotated relative to the world frame
angle = np.pi / 8
R2_world = np.array([
    [np.cos(angle), 0, np.sin(angle)],
    [0, 1, 0],
    [-np.sin(angle), 0, np.cos(angle)]
])
T2_world = np.array([-1., 0., 0.5])
C2_world = -R2_world.T @ T2_world # Camera center in world frame

# Common intrinsic matrix K
K = np.array([
    [1000, 0, 500],
    [0, 1000, 500],
    [0, 0, 1]
])
K_inv = np.linalg.inv(K)

# 3. Project the 3D point onto each camera's image plane
# Project to Camera 1
P_cam1 = R1_world @ P_world + T1_world
p_img1_homogeneous = K @ P_cam1
p_img1 = (p_img1_homogeneous / p_img1_homogeneous[2])[:2]

# Project to Camera 2
P_cam2 = R2_world @ P_world + T2_world
p_img2_homogeneous = K @ P_cam2
p_img2 = (p_img2_homogeneous / p_img2_homogeneous[2])[:2]

# 4. Form 3D lines (rays) from the image points
# For Camera 1, its frame is the world frame.
ray1_dir_cam1 = K_inv @ np.array([p_img1[0], p_img1[1], 1])
ray1_dir_world = R1_world.T @ ray1_dir_cam1 # In this case, R1 is identity
ray1_origin_world = C1_world

# For Camera 2, the ray is initially in Camera 2's local frame.
ray2_dir_cam2 = K_inv @ np.array([p_img2[0], p_img2[1], 1])
ray2_origin_cam2 = np.array([0., 0., 0.]) # Origin of camera 2's frame

# 5. The CRUCIAL TRANSFORMATION STEP
# We MUST transform Camera 2's ray into the world frame to match Camera 1's ray.
# Plücker coordinates or any other line representation would require this same step.
print("Demonstrating the necessary transformation:")
print("Ray direction from Cam 2 (in its own local frame):", np.round(ray2_dir_cam2, 3))
ray2_dir_world = R2_world.T @ ray2_dir_cam2
print("Ray direction from Cam 2 (transformed to world frame):", np.round(ray2_dir_world, 3))
print("-" * 30)

# The origin of the ray from Camera 2 is its center, which is already in world coordinates.
ray2_origin_world = C2_world

# 6. Triangulate the point using the two rays now in the SAME (world) frame
reconstructed_P = triangulate_point(
    ray1_origin_world, ray1_dir_world,
    ray2_origin_world, ray2_dir_world
)

print("Original 3D point:        ", P_world)
print("Reconstructed 3D point:   ", np.round(reconstructed_P, 3))
print("\nThe successful reconstruction was only possible after transforming the ray from")
print("the second camera's reference frame into the common world reference frame.")
print("This shows that Plücker coordinates cannot directly yield a solution; a transformation is needed.")
