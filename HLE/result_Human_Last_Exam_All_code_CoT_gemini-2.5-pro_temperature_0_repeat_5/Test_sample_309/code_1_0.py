def calculate_max_marginalized_landmarks(N, M):
    """
    Calculates the maximum number of landmarks that can be marginalized in a
    bundle adjustment problem while ensuring the system remains solvable.

    Args:
        N (int): The total number of landmarks.
        M (int): The total number of cameras (for context, not used in the final formula).
    """
    print(f"Given a system with {N} landmarks and {M} cameras.")
    print("To ensure the bundle adjustment problem is solvable, we must fix the 7-degree-of-freedom gauge ambiguity (3 translation, 3 rotation, 1 scale).")
    
    # A minimal set of 3D points is required to define a stable 3D coordinate frame.
    # 1st point fixes the origin (3 DOF).
    # 2nd point fixes an axis and scale (3 DOF).
    # 3rd point fixes the final rotation (1 DOF).
    min_landmarks_to_keep = 3
    
    print(f"A minimum of {min_landmarks_to_keep} landmarks must be kept in the optimization state to serve as an anchor for the coordinate system.")
    print("These anchor landmarks cannot be marginalized.")
    
    # The maximum number of landmarks that can be marginalized is the total
    # number of landmarks minus the minimum number that must be kept.
    max_marginalized = N - min_landmarks_to_keep
    
    print("\nThe calculation for the maximum number of marginalized landmarks is:")
    # The final output prints the equation with the specific numbers
    print(f"Max Marginalized = (Total Landmarks) - (Minimum Landmarks to Keep)")
    print(f"Max Marginalized = {N} - {min_landmarks_to_keep} = {max_marginalized}")
    print("\nTherefore, the general formula is N - 3.")

# --- Main execution ---
# We can use any example values for N and M. The formula N - 3 is independent of the specific values.
if __name__ == "__main__":
    # Example values from a hypothetical scenario
    total_landmarks = 1000
    total_cameras = 50
    calculate_max_marginalized_landmarks(total_landmarks, total_cameras)