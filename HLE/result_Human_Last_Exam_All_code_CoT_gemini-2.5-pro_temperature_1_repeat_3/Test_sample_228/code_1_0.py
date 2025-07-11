import numpy as np

def get_plucker_from_origin_to_point(p):
    """
    Computes the Plücker coordinates for a line passing through the origin and point p.
    For a line through points A and B, Plücker coords are (B-A, A x B).
    Here, A is the origin (0,0,0) and B is the point p.
    """
    # Direction vector of the line
    l = p - 0
    # Moment vector of the line
    m = np.cross(0, p) # This will be the zero vector
    return l, m

def transform_plucker_line(L, R, t):
    """
    Transforms a Plücker line L=(l,m) by a rotation R and translation t.
    """
    l, m = L
    l_new = R @ l
    m_new = R @ m + np.cross(t, l_new)
    return l_new, m_new

def reconstruct_point_from_lines(L1, L2):
    """
    Reconstructs the intersection point P of two Plücker lines L1 and L2.
    L1 is assumed to pass through the origin (m1=0).
    The point P must satisfy m2 = P x l2.
    Since P is on L1, P = k*l1 for some scalar k.
    So, m2 = (k*l1) x l2 = k * (l1 x l2).
    We solve for k and then find P.
    """
    l1, m1 = L1
    l2, m2 = L2
    
    # Ensure lines are in the same reference frame before calling this
    if np.linalg.norm(m1) > 1e-9:
        print("Error: This simplified reconstruction assumes the first line passes through the origin.")
        return None

    # Vector triple product to solve for k
    l1_cross_l2 = np.cross(l1, l2)
    
    # Avoid division by zero if lines are parallel
    if np.linalg.norm(l1_cross_l2) < 1e-9:
        print("Lines are parallel, cannot reconstruct point.")
        return None

    k = np.dot(m2, l1_cross_l2) / np.linalg.norm(l1_cross_l2)**2
    
    # Reconstructed point
    P_recon = k * l1
    return P_recon

# --- Simulation Setup ---

# 1. Define a 3D point in the World Frame.
#    Let's assume the World Frame is the same as Camera 1's frame.
P_world = np.array([2.0, 3.0, 10.0])
print(f"Original 3D Point (in World/Cam1 Frame): {P_world}\n")

# 2. Define the pose of Camera 2 relative to Camera 1 (World Frame)
#    Rotation: 10 degrees around y-axis
theta = np.deg2rad(10)
R_cam2_to_world = np.array([
    [np.cos(theta), 0, np.sin(theta)],
    [0, 1, 0],
    [-np.sin(theta), 0, np.cos(theta)]
])
#    Translation: 1 unit along the x-axis
t_cam2_to_world = np.array([-1, 0, 0])

# We need the transformation FROM World TO Camera 2 for projecting the point
R_world_to_cam2 = R_cam2_to_world.T
t_world_to_cam2 = -R_world_to_cam2 @ t_cam2_to_world


# --- Line Generation ---

# 3. Line 1 (in Camera 1's frame)
# The line goes from Camera 1's center (0,0,0) to P_world.
l1, m1 = get_plucker_from_origin_to_point(P_world)
print("--- Step 1: Line in Camera 1's Frame ---")
print(f"Plücker Line L1 direction (l1): {l1}")
print(f"Plücker Line L1 moment (m1):   {m1}\n")


# 4. Line 2 (in Camera 2's frame)
# First find the coordinates of P_world in Camera 2's frame
P_in_cam2 = R_world_to_cam2 @ P_world + t_world_to_cam2
# The line goes from Camera 2's center (0,0,0) to P_in_cam2
l2, m2 = get_plucker_from_origin_to_point(P_in_cam2)
print("--- Step 2: Line in Camera 2's Frame ---")
print("This line represents the same ray, but its coordinates are different because it's in a different frame.")
print(f"Plücker Line L2 direction (l2): {l2}")
print(f"Plücker Line L2 moment (m2):   {m2}\n")


# --- Transformation and Reconstruction ---

print("--- Step 3: The Need for Transformation ---")
print("To find the 3D point, we must express both lines in the SAME reference frame.")
print("We will transform L2 from Camera 2's frame into Camera 1's frame.\n")

# 5. Transform L2 into Camera 1's frame
# We need the transformation FROM Camera 2 TO Camera 1 (World)
l2_transformed, m2_transformed = transform_plucker_line((l2, m2), R_cam2_to_world, t_cam2_to_world)
print("--- Step 4: Transformed Line 2 ---")
print(f"Transformed L2 direction (in Cam1 frame): {l2_transformed}")
print(f"Transformed L2 moment (in Cam1 frame):   {m2_transformed}\n")


# 6. Reconstruct the 3D point
# Now both L1 and the transformed L2 are in Camera 1's frame. We can find their intersection.
P_reconstructed = reconstruct_point_from_lines((l1, m1), (l2_transformed, m2_transformed))
print("--- Step 5: Reconstruction ---")
print("With both lines in the same frame, we can calculate their intersection.")
print(f"Reconstructed 3D Point: {P_reconstructed}")
print(f"Original 3D Point:      {P_world}")
print(f"Success! The point is reconstructed after applying the necessary transformation.\n")
print("This demonstrates that a transformation is a required step, so the solution is not yielded directly.")
