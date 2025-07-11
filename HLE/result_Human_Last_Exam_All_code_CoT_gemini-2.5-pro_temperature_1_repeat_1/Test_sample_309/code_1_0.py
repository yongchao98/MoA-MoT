def calculate_max_marginalized_landmarks(N, M):
    """
    This function calculates the maximum number of landmarks that can be
    marginalized in a bundle adjustment problem while ensuring solvability.

    Args:
      N (int): The total number of landmarks in the problem.
      M (int): The total number of cameras in the problem. (Note: M is not
               needed for this specific calculation but is part of the problem's context).
    """

    print("Analyzing the Bundle Adjustment (BA) problem:")
    print(f"Given: {N} landmarks and {M} cameras.")
    print("Strategy: Marginalize landmarks to solve for cameras first using the Schur complement.")
    print("-" * 50)
    print("Key Insight: The BA problem has a 7-DOF gauge ambiguity (3 translation, 3 rotation, 1 scale).")
    print("To make the problem solvable, we must fix this gauge by anchoring the coordinate system.")
    print("\nA common and robust method is to anchor the world frame using a minimal set of landmarks.")
    
    # The number of landmarks required to define a 3D coordinate frame.
    anchor_landmarks = 3
    
    print(f"A minimum of {anchor_landmarks} landmarks are required to define a 3D frame and fix the gauge:")
    print("  - Landmark 1 fixes the origin (3 DOFs).")
    print("  - Landmark 2 fixes an axis (2 DOFs).")
    print("  - Landmark 3 fixes a plane (1 DOF).")
    print("  - One distance is fixed to set the scale (1 DOF).")
    print("\nThese anchor landmarks cannot be marginalized because they are needed to apply the constraints.")
    print("-" * 50)
    
    # The maximum number of landmarks that can be marginalized is the total
    # number of landmarks minus the number of anchor landmarks.
    max_marginalized_landmarks = N - anchor_landmarks
    
    print("Calculation:")
    print(f"Maximum number of marginalized landmarks = (Total Landmarks) - (Anchor Landmarks)")
    print(f"Result = {N} - {anchor_landmarks}")
    print(f"Result = {max_marginalized_landmarks}")


# Let's use example values for N and M.
# The question is symbolic, so the exact values are not important for the formula.
# We assume N is large enough (N >= 3).
N_landmarks = 100
M_cameras = 20

calculate_max_marginalized_landmarks(N_landmarks, M_cameras)