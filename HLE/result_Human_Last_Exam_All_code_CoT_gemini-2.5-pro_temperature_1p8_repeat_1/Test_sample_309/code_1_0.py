import argparse

def solve_bundle_adjustment_marginalization(N, M):
    """
    Calculates the maximum number of landmarks that can be marginalized in a
    bundle adjustment problem while ensuring the problem remains solvable.

    Args:
        N (int): The total number of landmarks.
        M (int): The total number of cameras.
    """

    print("Step 1: Understanding the problem of Bundle Adjustment (BA).")
    print("BA optimizes both 3D landmark positions (structure) and camera poses (motion).")
    print(f"We have N = {N} landmarks and M = {M} cameras.\n")

    print("Step 2: Identifying the challenge - Gauge Ambiguity.")
    print("A BA system has a 7-degree-of-freedom (DOF) ambiguity. We can move, rotate,")
    print("and scale the entire scene without changing the measurements (image projections).")
    print("To find a unique solution, we must 'anchor' or 'fix' the coordinate system.\n")
    
    print("Step 3: Anchoring the system using landmarks.")
    print("Fixing a camera's pose resolves 6 DOF (3 rotation, 3 translation).")
    print("The remaining 1 DOF (scale) and other rotational ambiguities must be resolved by")
    print("keeping a minimum number of landmarks in the optimization.\n")
    
    print("   - Keeping 1 landmark is not enough to fix scale.")
    print("   - Keeping 2 landmarks can fix scale, but allows rotation around the axis connecting them.")
    print("   - Keeping 3 non-collinear landmarks defines a plane, which fully anchors the")
    print("     geometry and resolves all ambiguities.\n")

    # The minimum number of landmarks that must be kept in the state vector.
    min_landmarks_to_keep = 3
    print(f"Step 4: Calculating the maximum number of landmarks to marginalize.")
    print(f"To ensure the problem remains solvable, a minimum of {min_landmarks_to_keep} landmarks must be kept.")
    
    if N < min_landmarks_to_keep:
        print("\nError: The total number of landmarks is less than the minimum required to solve the problem.")
        return

    # Calculate the maximum number of landmarks that can be marginalized.
    max_marginalized_landmarks = N - min_landmarks_to_keep

    print("\nThe maximum number of landmarks that can be marginalized is the total number")
    print("of landmarks minus the minimum number we must keep.\n")
    
    # Final equation output, as requested.
    print("--- Final Calculation ---")
    print(f"Maximum Marginalized Landmarks = (Total Landmarks) - (Minimum Landmarks to Keep)")
    print(f"Maximum Marginalized Landmarks = {N} - {min_landmarks_to_keep} = {max_marginalized_landmarks}")
    print("-------------------------")


if __name__ == '__main__':
    # Example values for N and M.
    # The value of M is not directly used in the final calculation but is part of the problem statement.
    example_N = 100
    example_M = 10
    
    solve_bundle_adjustment_marginalization(example_N, example_M)