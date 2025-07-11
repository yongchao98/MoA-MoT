import numpy as np

def skew_symmetric_matrix(v):
    """Creates a skew-symmetric matrix from a 3-vector."""
    return np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

def triangulate_midpoint(p1, d1, p2, d2):
    """
    Finds the midpoint of the shortest segment between two skew lines.
    Line 1 is defined by point p1 and direction d1.
    Line 2 is defined by point p2 and direction d2.
    """
    # Direction of the common perpendicular
    d_perp = np.cross(d1, d2)
    d_perp_norm = np.linalg.norm(d_perp)
    if d_perp_norm < 1e-6:
        # Lines are parallel, cannot find a unique midpoint this way.
        # Returning the midpoint between the initial points as a fallback.
        return (p1 + p2) / 2

    # Solve for the parameters t1 and t2 where the lines are closest
    A = np.array([d1, -d2, d_perp]).T
    b = p2 - p1
    try:
        # We only need the first two parameters for the points on the lines
        params = np.linalg.solve(A, b)
        t1 = params[0]
        t2 = params[1]

        # Points on each line that are closest to the other line
        closest_p1 = p1 + t1 * d1
        closest_p2 = p2 + t2 * d2

        # The triangulated point is the midpoint of the segment connecting them
        triangulated_point = (closest_p1 + closest_p2) / 2
        return triangulated_point
    except np.linalg.LinAlgError:
        # Matrix is singular (e.g., lines are parallel or intersecting)
        return (p1 + p2) / 2 # Fallback


# 1. Setup the Scene
# The world coordinate system is the same as Camera 1's coordinate system.
P_world = np.array([2.0, 3.0, 10.0])

# Camera 1 (at the origin, looking along Z-axis)
C1 = np.array([0.0, 0.0, 0.0]) # Center
R1 = np.identity(3)
t1 = -R1 @ C1
P1 = np.hstack([R1, t1.reshape(-1, 1)]) # Extrinsic matrix

# Camera 2 (stereo setup, translated along X-axis)
baseline = 5.0
C2 = np.array([baseline, 0.0, 0.0]) # Center
R2 = np.identity(3)
t2 = -R2 @ C2
P2 = np.hstack([R2, t2.reshape(-1, 1)]) # Extrinsic matrix

# 2. Project and Add Noise
# Project point onto normalized image planes (K=Identity)
P_homogeneous = np.append(P_world, 1)
p1_proj = P1 @ P_homogeneous
p1_img = p1_proj[:2] / p1_proj[2]

p2_proj = P2 @ P_homogeneous
p2_img = p2_proj[:2] / p2_proj[2]

# Add a small amount of noise to simulate real measurements
noise_level = 0.001
p1_noisy = p1_img + np.random.normal(0, noise_level, 2)
p2_noisy = p2_img + np.random.normal(0, noise_level, 2)

# 3. Back-project to Rays (in the World/Camera 1 frame)
# Ray 1
# Direction vector from camera center to the point on the image plane
d1 = np.array([p1_noisy[0], p1_noisy[1], 1.0])
d1 = d1 / np.linalg.norm(d1) # Normalize direction

# Ray 2
# Back-project point in Cam2's frame, then transform direction to world frame
d2_cam2 = np.array([p2_noisy[0], p2_noisy[1], 1.0])
d2 = R2.T @ d2_cam2 # R2 is identity, so this doesn't change anything
d2 = d2 / np.linalg.norm(d2) # Normalize direction

# 4. Create Plücker Coordinates for the rays
# L = (direction, moment) where moment = point x direction
m1 = np.cross(C1, d1)
L1 = (d1, m1)

m2 = np.cross(C2, d2)
L2 = (d2, m2)

# 5. Check for Intersection using the Reciprocal Product
# The product is d1.m2 + d2.m1
reciprocal_product = np.dot(d1, m2) + np.dot(d2, m1)

# 6. Triangulate
# Since they don't intersect, find the point of closest approach
P_triangulated = triangulate_midpoint(C1, d1, C2, d2)

# 7. Print Results
print("This script demonstrates 3D triangulation and its limitations.")
print(f"The original 3D point is at: [{P_world[0]:.4f}, {P_world[1]:.4f}, {P_world[2]:.4f}]")
print("-" * 60)
print(f"After adding noise, the two back-projected rays are skew (non-intersecting).")
print(f"We can prove this with the Plucker reciprocal product.")
print(f"If the lines intersect, the product is 0.")
print(f"Calculated reciprocal product: {reciprocal_product:.8f}")
print("-" * 60)
print("Because the lines do not intersect, a 'direct' triangulation is not possible.")
print("Instead, we find an optimal point, such as the midpoint of the segment of")
print("the shortest distance between the two lines.")
print(f"The final triangulated point is: [{P_triangulated[0]:.4f}, {P_triangulated[1]:.4f}, {P_triangulated[2]:.4f}]")
print("(Note: The result is in the reference frame of Camera 1, which we defined as the world frame)")
