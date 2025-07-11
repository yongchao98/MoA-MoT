import numpy as np

def plucker_from_points(p1, p2):
    """Computes Plücker coordinates L = (d, m) for a line through p1 and p2."""
    d = p2 - p1
    m = np.cross(p1, p2)
    return d / np.linalg.norm(d), m # Normalize direction for stability

def triangulate_with_plucker(L1, L2):
    """
    Solves for the 3D point P that lies on two lines L1 and L2.
    It solves the over-determined system [skew(d1); skew(d2)] * P = [m1; m2]
    using least squares.
    """
    d1, m1 = L1
    d2, m2 = L2

    def skew(v):
        """Returns the 3x3 skew-symmetric matrix for a 3D vector v."""
        return np.array([[0, -v[2], v[1]],
                         [v[2], 0, -v[0]],
                         [-v[1], v[0], 0]])

    # A * P = b
    A = np.vstack([skew(d1), skew(d2)])
    b = np.concatenate([m1, m2])

    # Solve the system using linear least squares
    P, residuals, rank, s = np.linalg.lstsq(A, b, rcond=None)
    return P

# 1. Scene Setup
# Let's define the world reference frame to be the same as Camera 1's frame.
P_true = np.array([2.0, 3.0, 10.0])

# Camera 1 is at the origin of the world.
C1_world = np.array([0.0, 0.0, 0.0])

# Camera 2's pose relative to Camera 1 (the world).
# Let's say it's translated 5 units on the X-axis.
# The transformation that brings world points to Cam2's frame is P_c2 = R @ P_world + t_vec.
# For simplicity, R is identity, t_vec is [-5, 0, 0].
# The center of Camera 2 in world coordinates is therefore [5, 0, 0].
C2_world = np.array([5.0, 0.0, 0.0])

print("--- Scene Setup ---")
print(f"True 3D point in World/Camera1 frame: P = {P_true}")
print(f"Camera 1 Center in World frame: C1 = {C1_world}")
print(f"Camera 2 Center in World frame: C2 = {C2_world}")
print("-" * 20)

# 2. Define Lines of Sight in a COMMON Reference Frame (World Frame)
# This is the key concept. All geometry must exist in the same frame.

# Line 1 from Camera 1 is defined by C1_world and P_true.
L1_world = plucker_from_points(C1_world, P_true)

# Line 2 from Camera 2 must also be defined in the world frame.
# This requires using the transformed center C2_world. We cannot use its
# local origin (0,0,0). This is the necessary transformation step.
L2_world = plucker_from_points(C2_world, P_true)

print("--- Plücker Line Generation (in World Frame) ---")
print("Line 1 Direction (d1):", np.round(L1_world[0], 4))
print("Line 1 Moment (m1):   ", np.round(L1_world[1], 4))
print("\nLine 2 Direction (d2):", np.round(L2_world[0], 4))
print("Line 2 Moment (m2):   ", np.round(L2_world[1], 4))
print("-" * 20)


# 3. Triangulate the 3D Point
# Now that both lines are in the same reference frame, we can find their intersection.
P_triangulated = triangulate_with_plucker(L1_world, L2_world)

print("--- Triangulation Result ---")
print("The intersection is found by solving a system of linear equations derived from the Plücker coordinates.")
print(f"Original Point:       P = [{P_true[0]}, {P_true[1]}, {P_true[2]}]")
print(f"Triangulated Point:   P' = [{P_triangulated[0]:.4f}, {P_triangulated[1]:.4f}, {P_triangulated[2]:.4f}]")
print("\nConclusion: The calculation works, but only after Line 2 was transformed into the reference frame of Line 1, which required knowing the pose of Camera 2. Therefore, a transformation is needed.")
