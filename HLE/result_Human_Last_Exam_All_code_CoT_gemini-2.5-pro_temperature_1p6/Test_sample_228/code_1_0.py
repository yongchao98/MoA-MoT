import numpy as np

def skew_symmetric_matrix(v):
    """Creates a skew-symmetric matrix from a 3D vector."""
    return np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

def get_plucker_line(p1, p2):
    """
    Calculates the Plücker coordinates for a line passing through points p1 and p2.
    The points must be in the same reference frame.
    """
    d = p2 - p1
    d = d / np.linalg.norm(d)  # Normalize the direction vector
    m = np.cross(p1, d)
    return d, m

def main():
    """
    Demonstrates triangulation of a 3D point using Plücker coordinates,
    with all calculations performed directly in one camera's reference frame.
    """
    print("--- 3D Triangulation using Plücker Coordinates in a Camera's Reference Frame ---")

    # Step 1: Define the scene in Camera 1's reference frame.
    # Camera 1 is at the origin (0,0,0) of this frame.
    cam1_center = np.array([0.0, 0.0, 0.0])

    # The ground truth 3D point we want to find.
    # Its coordinates are given in Camera 1's frame.
    X_ground_truth = np.array([1.0, 2.0, 5.0])
    print(f"\nGround Truth 3D Point (in Cam 1's frame): {X_ground_truth}")

    # Define the pose of Camera 2 relative to Camera 1.
    # This transformation takes points from Cam 1's frame to Cam 2's frame.
    R = np.array([[0.995, 0.0, 0.0998], [0.0, 1.0, 0.0], [-0.0998, 0.0, 0.995]]) # Rotation
    t = np.array([-2.0, 0.0, 0.0]) # Translation

    # Calculate the center of Camera 2 in Camera 1's frame.
    # p_c2 = R @ p_c1 + t => p_c1 = inv(R) @ (p_c2 - t)
    # If p_c2 is the origin of cam2 (0,0,0), then p_c1 is the position in cam1's frame.
    cam2_center = -np.linalg.inv(R) @ t
    print(f"Camera 2 Center (in Cam 1's frame): {cam2_center}\n")

    # Step 2: Define the two viewing rays (lines) in Camera 1's frame.
    # A viewing ray is a line from a camera center to the 3D point.
    print("Calculating Plücker coordinates for both viewing rays in Camera 1's frame...")
    # Ray 1: From Camera 1 center to the 3D point
    d1, m1 = get_plucker_line(cam1_center, X_ground_truth)
    print(f"  Line 1 Direction (d1): [{d1[0]:.4f} {d1[1]:.4f} {d1[2]:.4f}]")
    print(f"  Line 1 Moment   (m1): [{m1[0]:.4f} {m1[1]:.4f} {m1[2]:.4f}]")

    # Ray 2: From Camera 2 center to the 3D point
    d2, m2 = get_plucker_line(cam2_center, X_ground_truth)
    print(f"  Line 2 Direction (d2): [{d2[0]:.4f} {d2[1]:.4f} {d2[2]:.4f}]")
    print(f"  Line 2 Moment   (m2): [{m2[0]:.4f} {m2[1]:.4f} {m2[2]:.4f}]")

    # Step 3: Triangulate the point using the Plücker coordinates.
    # A point X on a line L=(d,m) satisfies the equation: d x X = m
    # This can be written as [d]_x X = m, where [d]_x is the skew-symmetric matrix of d.
    # We have two such equations, one for each line, forming an overdetermined system.
    # A = [[d1]_x]   b = [m1]
    #     [[d2]_x]       [m2]
    # We solve the linear system A*X = b for X using a least-squares method.

    A1 = skew_symmetric_matrix(d1)
    A2 = skew_symmetric_matrix(d2)
    A = np.vstack([A1, A2])
    b = np.concatenate([m1, m2])

    # Solve for X using the pseudo-inverse: X = (A^T A)^-1 A^T b
    X_triangulated = np.linalg.pinv(A) @ b

    # Step 4: Output the final result.
    # The result is directly in Camera 1's coordinate frame because all inputs were in that frame.
    print("\n--- Triangulation Result ---")
    print("The final equation solved is A * X = b, where:")
    print(f"A (6x3 matrix from d1, d2):\n{A}")
    print(f"b (6x1 vector from m1, m2):\n{b.reshape(-1, 1)}")
    print(f"\nFinal Triangulated Point (in Cam 1's frame): [{X_triangulated[0]:.4f} {X_triangulated[1]:.4f} {X_triangulated[2]:.4f}]")

if __name__ == '__main__':
    main()