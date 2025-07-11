def solve_bundle_adjustment_marginalization():
    """
    This function determines the maximum number of landmarks that can be marginalized
    in a bundle adjustment problem while ensuring the system remains solvable.

    The reasoning is as follows:
    1. A standard monocular bundle adjustment problem has a 7 degree-of-freedom (DoF)
       gauge ambiguity (3 for translation, 3 for rotation, 1 for scale).
    2. To make the problem solvable, this ambiguity must be removed by adding at least
       7 constraints, typically by fixing a reference frame.
    3. If we use landmarks to define this frame, we must determine the minimum number
       of landmarks required to provide these 7 constraints.
    4. A minimum of 3 landmarks are needed:
       - The 1st landmark can fix translation (3 DoF).
       - The 2nd and 3rd landmarks can fix rotation (3 DoF).
       - A specified distance involving at least two of these landmarks can fix scale (1 DoF).
    5. These 3 landmarks, which are used to anchor the solution, cannot be marginalized.
       Their parameters must remain in the optimization to be constrained.
    6. Given N total landmarks, the rest of the landmarks can be marginalized.
    """

    # N represents the total number of landmarks. It is treated as a symbol.
    n_landmarks_symbol = "N"

    # This is the minimum number of landmarks that must be kept in the optimization
    # to fix the 7-DoF gauge freedom.
    min_landmarks_to_keep = 3

    # The final equation for the maximum number of marginalized landmarks is N - 3.
    # The instructions require outputting each number in the final equation.
    # Here, the number is 3.
    print(f"The maximum number of landmarks that can be marginalized is: {n_landmarks_symbol} - {min_landmarks_to_keep}")

# Execute the function to print the result.
solve_bundle_adjustment_marginalization()