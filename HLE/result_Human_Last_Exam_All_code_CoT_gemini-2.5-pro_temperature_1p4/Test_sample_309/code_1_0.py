import argparse

def solve(N, M):
    """
    Calculates the maximum number of landmarks that can be marginalized in a bundle adjustment problem.

    Args:
        N (int): The total number of landmarks.
        M (int): The total number of cameras.
    """
    # The problem of solving for camera poses in a bundle adjustment setting
    # is only well-posed if the underlying 3D structure is rigid and defines a
    # unique coordinate frame (up to a global similarity transform).
    # Geometrically, a minimum of 3 non-collinear landmarks are required to
    # create this stable 3D reference frame.

    # While it's algebraically possible to marginalize all N landmarks (if N >= 3),
    # a common and robust strategy is to keep a minimal set of landmarks as "anchors"
    # in the main optimization problem to ensure the system remains well-conditioned
    # and solvable. These anchor landmarks are not marginalized.

    # The minimum number of landmarks required to form this anchor frame is 3.
    num_anchor_landmarks = 3

    # The problem statement implies we must have a solvable system, which means
    # N must be at least 3. If N is less than 3, the problem of solving for
    # camera geometry is ill-posed.
    if N < num_anchor_landmarks:
        print(f"With only {N} landmarks, the problem is not geometrically solvable for camera poses.")
        # In this case, one cannot even form a basis, so 0 landmarks can be marginalized.
        # The result of N-3 would be negative, which doesn't make sense.
        # The most logical answer in this case is 0.
        max_marginalized = 0
    else:
        # If we have N total landmarks, and we must designate 3 as non-marginalized anchors
        # to ensure the geometric stability of the camera optimization problem, then the
        # remaining landmarks can be marginalized out using the Schur complement.
        max_marginalized = N - num_anchor_landmarks

    print(f"Given {N} landmarks and {M} cameras:")
    print("To ensure the problem remains geometrically solvable, a stable reference frame is needed.")
    print(f"A minimum of {num_anchor_landmarks} landmarks are required to form this reference frame.")
    print("These anchor landmarks are kept in the optimization state and are not marginalized.")
    print(f"The remaining {N} - {num_anchor_landmarks} landmarks can be marginalized.")
    print("\nFinal Equation:")
    print(f"Maximum Marginalized Landmarks = Total Landmarks (N) - Anchor Landmarks")
    print(f"Maximum Marginalized Landmarks = {N} - {num_anchor_landmarks} = {max_marginalized}")

if __name__ == '__main__':
    # The problem is a theoretical one, so N and M are symbolic.
    # We can use example values to run the code.
    # Let's use a hypothetical N=50 landmarks and M=10 cameras.
    parser = argparse.ArgumentParser(description="Calculate max marginalized landmarks in BA.")
    parser.add_argument('--N', type=int, default=50, help='Total number of landmarks')
    parser.add_argument('--M', type=int, default=10, help='Total number of cameras')
    args = parser.parse_args()
    solve(args.N, args.M)
