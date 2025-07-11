import sys

def calculate_max_marginalized_landmarks(N, M):
    """
    Calculates the maximum number of landmarks that can be marginalized in a
    bundle adjustment problem while keeping it solvable.

    Args:
        N (int): The total number of landmarks.
        M (int): The total number of cameras.
    """
    print(f"Given the problem with {N} landmarks and {M} cameras.")
    print("To ensure the bundle adjustment problem remains solvable, we need a stable geometric reference frame.")
    print("A minimum of 3 non-collinear landmarks is required to define a stable 3D reference frame and resolve geometric ambiguities.")
    print("\nLet 'k' be the number of landmarks that are marginalized.")
    print("The number of landmarks remaining in the active optimization state is (N - k).")
    print("To maintain a solvable system, this number must be at least 3.")
    print("\nThis gives us the inequality:")
    print("N - k >= 3")
    print("\nSolving for k, we find the maximum number of landmarks that can be marginalized:")
    print("k <= N - 3")

    if N < 3:
        max_k = 0
        print("\nSince N is less than 3, the problem is not well-posed to begin with.")
        print("Result: No landmarks can be marginalized.")
        print("Maximum number of marginalized landmarks = 0")
    else:
        max_k = N - 3
        print("\nFinal Calculation:")
        print(f"Maximum marginalized landmarks = N - 3 = {N} - 3 = {max_k}")


if __name__ == '__main__':
    # You can change these values to test with different numbers of landmarks and cameras.
    # We get these from command line arguments if provided, otherwise use defaults.
    try:
        num_landmarks = int(sys.argv[1]) if len(sys.argv) > 1 else 50
        num_cameras = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    except (ValueError, IndexError):
        print("Usage: python your_script_name.py [num_landmarks] [num_cameras]")
        print("Using default values: N=50, M=10")
        num_landmarks = 50
        num_cameras = 10

    calculate_max_marginalized_landmarks(num_landmarks, num_cameras)