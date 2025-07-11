import sys

def solve_bundle_adjustment_marginalization(N, M):
    """
    Calculates the maximum number of landmarks that can be marginalized in a
    bundle adjustment problem under the assumption that 3 landmarks are needed
    to fix the gauge freedom.

    Args:
        N (int): The number of landmarks.
        M (int): The number of cameras.
    """
    if not isinstance(N, int) or not isinstance(M, int):
        print("Error: N (number of landmarks) and M (number of cameras) must be integers.")
        return

    if N < 3:
        print(f"Error: With N={N} landmarks, you cannot keep 3 landmarks to establish a coordinate system.")
        print("The problem is not solvable under this assumption.")
        return
        
    # To ensure the problem remains solvable, we need to fix the 7-DOF gauge freedom.
    # A standard way to do this is to anchor the coordinate system to a set of static landmarks.
    # A minimum of 3 non-collinear landmarks are required to robustly define a 3D coordinate frame
    # (origin, orientation, and scale).
    # These 3 landmarks must remain in the active set of variables and cannot be marginalized.
    # Therefore, the maximum number of landmarks that can be marginalized is N - 3.
    
    landmarks_to_keep = 3
    max_marginalized_landmarks = N - landmarks_to_keep

    print(f"Given:")
    print(f"  Number of landmarks (N) = {N}")
    print(f"  Number of cameras (M)   = {M}")
    print("\nTo ensure the problem is solvable, we must fix the 7-DOF gauge freedom (translation, rotation, scale).")
    print("A common method is to define the world coordinate frame using a set of anchor landmarks.")
    print(f"A minimum of {landmarks_to_keep} non-collinear landmarks are required to define a stable 3D frame.")
    print("These landmarks must not be marginalized.")
    print("\nFinal Calculation:")
    print(f"Maximum number of marginalized landmarks = (Total Landmarks) - (Landmarks to Keep)")
    print(f"Result = {N} - {landmarks_to_keep} = {max_marginalized_landmarks}")


# --- Example Usage ---
# You can change these values to test with different numbers
# Default values for the example
num_landmarks = 100
num_cameras = 10

# You can also provide command-line arguments for N and M
if len(sys.argv) == 3:
    try:
        num_landmarks = int(sys.argv[1])
        num_cameras = int(sys.argv[2])
    except ValueError:
        print("Invalid command-line arguments. Please provide two integers for N and M.")
        sys.exit(1)

solve_bundle_adjustment_marginalization(num_landmarks, num_cameras)