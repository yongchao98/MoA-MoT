import numpy as np

def run_triangulation_demo():
    """
    Demonstrates that 3D triangulation yields a point in the world reference frame,
    and a transformation is required to express it in a camera's local frame.
    """

    # --- 1. Define a World Frame and a 3D Point ---
    # The world frame is our global coordinate system.
    # We place a 3D point 'P' in this world frame.
    P_world = np.array([4.0, 2.0, 15.0])
    print(f"Original 3D point in world frame (P_world): {P_world}\n")

    # --- 2. Define Two Cameras in the World Frame ---
    # The position and orientation of each camera are its 'extrinsic' parameters.

    # Camera 1 is at the world origin.
    C1_world = np.array([0.0, 0.0, 0.0])
    R1 = np.identity(3)  # No rotation relative to the world

    # Camera 2 is shifted along the X-axis.
    C2_world = np.array([10.0, 0.0, 0.0])
    R2 = np.identity(3)  # No rotation relative to the world

    print(f"Camera 1 position in world frame (C1_world): {C1_world}")
    print(f"Camera 2 position in world frame (C2_world): {C2_world}\n")

    # --- 3. Define Viewing Rays in the World Frame ---
    # A viewing ray is a line from a camera center to the 3D point.
    # The direction vectors are calculated in the world frame.
    # While Plücker coordinates are a valid representation, for this demonstration,
    # using simple direction vectors is clearer and shows the same principle.
    
    # Normalize the direction vectors
    d1_world = (P_world - C1_world) / np.linalg.norm(P_world - C1_world)
    d2_world = (P_world - C2_world) / np.linalg.norm(P_world - C2_world)

    print(f"Ray 1 direction in world frame: {np.round(d1_world, 3)}")
    print(f"Ray 2 direction in world frame: {np.round(d2_world, 3)}\n")

    # --- 4. Triangulate the Point ---
    # We find the intersection of the two rays: C1 + a*d1 = C2 + b*d2
    # This solves for 'a' and 'b', which are the distances along each ray.
    # The calculation uses world-frame coordinates (C1, C2, d1, d2).
    A = np.vstack((d1_world, -d2_world)).T
    b = C2_world - C1_world
    
    # Solve the system of linear equations for [a, b]
    try:
        solution = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        # If the system is overdetermined (more equations than unknowns), use least squares
        solution = np.linalg.lstsq(A, b, rcond=None)[0]

    a = solution[0]

    # Calculate the triangulated point using the result.
    # The result MUST be in the same frame as the inputs: the world frame.
    P_triangulated = C1_world + a * d1_world
    print(f"Triangulated point coordinates: {np.round(P_triangulated, 3)}")
    print(f"Note: This is identical to the original P_world.\n")


    # --- 5. Find the Point's Coordinates in a Camera's Local Frame ---
    # To get the point's coordinates in a camera's frame, we must apply a transformation.
    # P_camera = R_transpose * (P_world - C_world)

    # Let's find the coordinates of the point relative to Camera 2.
    P_in_cam2_frame = R2.T @ (P_triangulated - C2_world)
    
    print("--- Final Analysis ---")
    print(f"The direct result of triangulation is the point P = {np.round(P_triangulated, 3)}")
    print("This result is in the WORLD reference frame.")
    print(f"To find the point's coordinates in Camera 2's local frame, a transformation is needed.")
    print(f"The transformed coordinates are: {np.round(P_in_cam2_frame, 3)}")
    print("\nAs you can see, the direct result of triangulation is not in the camera's reference frame.")
    print("Therefore, a different transformation is needed.")

if __name__ == '__main__':
    run_triangulation_demo()