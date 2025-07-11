import numpy as np

def skew_symmetric(v):
    """Creates a skew-symmetric matrix from a 3-vector."""
    return np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

def create_plucker_line(p1, p2):
    """Creates Plucker coordinates for a line passing through p1 and p2."""
    direction = p2 - p1
    direction /= np.linalg.norm(direction)
    moment = np.cross(p1, p2)
    # For a line passing through the origin (p1=0), moment is 0.
    # The standard representation is m = p1 x direction.
    moment = np.cross(p1, direction)
    return direction, moment

def transform_plucker_line(d_local, m_local, R, t):
    """Transforms a Plucker line from a local frame to a world frame.
    
    Args:
        d_local: Direction vector in local frame.
        m_local: Moment vector in local frame.
        R: Rotation matrix that transforms points from local to world.
        t: Translation vector (origin of local frame in world coords).
    """
    d_world = R @ d_local
    # The moment transformation formula is m_world = R @ m_local + t x (R @ d_local)
    m_world = R @ m_local + np.cross(t, d_world)
    return d_world, m_world

def triangulate_plucker(l1, l2):
    """Triangulates a 3D point from two Plucker lines."""
    d1, m1 = l1
    d2, m2 = l2
    
    # Build the system A*P = b for the overdetermined system
    # [d1]_x * P = m1
    # [d2]_x * P = m2
    A = np.vstack([skew_symmetric(d1), skew_symmetric(d2)])
    b = np.hstack([m1, m2])
    
    # Solve using the pseudo-inverse: P = (A^T * A)^-1 * A^T * b
    P, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    return P

# --- Simulation Setup ---

# 1. Define a 3D point in the world frame
P_world = np.array([2.0, 3.0, 10.0])

# 2. Define Camera 1 (at the world origin)
C1_pos = np.array([0.0, 0.0, 0.0]) # Position
C1_rot = np.identity(3)           # Orientation (identity)

# 3. Define Camera 2 (translated and rotated relative to world/Camera 1)
# Rotation around Y axis by -15 degrees
angle = np.deg2rad(-15)
C2_rot_world_to_cam = np.array([
    [np.cos(angle), 0, np.sin(angle)],
    [0, 1, 0],
    [-np.sin(angle), 0, np.cos(angle)]
])
# Transformation from camera frame to world frame is the inverse
C2_rot_cam_to_world = C2_rot_world_to_cam.T
C2_pos = np.array([4.0, 0.0, 0.0]) # Position of Cam2 in world frame

print(f"Original 3D Point P: {P_world}")
print("-" * 30)

# --- Line Generation and Transformation ---

# 4. Line 1 (from Camera 1's perspective)
# This line is already in the world frame because Camera 1 is at the origin.
print("Line 1 (Camera 1):")
d1_world, m1_world = create_plucker_line(C1_pos, P_world)
print(f"  - Plücker Coords (World Frame): d1={np.round(d1_world, 3)}, m1={np.round(m1_world, 3)}")

# 5. Line 2 (from Camera 2's perspective)
print("\nLine 2 (Camera 2):")
# First, find point P in Camera 2's local coordinate system
P_cam2 = C2_rot_world_to_cam @ (P_world - C2_pos)
# In Camera 2's frame, the line passes through its origin (0,0,0) and P_cam2
d2_cam2, m2_cam2 = create_plucker_line(np.array([0,0,0]), P_cam2)
print(f"  - Plücker Coords (Local Cam2 Frame): d2={np.round(d2_cam2, 3)}, m2={np.round(m2_cam2, 3)}")

# THIS IS THE CRUCIAL STEP: Transform Line 2 into the world frame
print("  - Applying transformation to bring Line 2 into World Frame...")
d2_world, m2_world = transform_plucker_line(d2_cam2, m2_cam2, C2_rot_cam_to_world, C2_pos)
print(f"  - Plücker Coords (World Frame): d2={np.round(d2_world, 3)}, m2={np.round(m2_world, 3)}")
print("-" * 30)

# --- Triangulation ---

# 6. Triangulate the point using the two lines now in the SAME (world) frame
L1_world = (d1_world, m1_world)
L2_world = (d2_world, m2_world)

P_reconstructed = triangulate_plucker(L1_world, L2_world)

print("\nTriangulation Result:")
print(f"The reconstructed 3D point is: {np.round(P_reconstructed, 3)}")
print("\nConclusion: The calculation required transforming the Plucker coordinates of the line from Camera 2's local frame into the world frame. Without this transformation, the triangulation would fail.")
print("Therefore, the solution is not yielded directly.")
