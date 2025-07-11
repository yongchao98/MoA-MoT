def solve_bundle_adjustment_question():
    """
    Explains the reasoning to find the maximum number of marginalizable landmarks in Bundle Adjustment.

    In a bundle adjustment problem with N landmarks and M cameras, to ensure a solvable system
    with a unique solution, we must first fix the inherent 7-degree-of-freedom gauge ambiguity.

    A standard method to do this is to fix the 3D coordinates of a minimal set of landmarks,
    which establishes a reference frame for the entire scene.

    The minimum number of non-collinear landmarks required to define a stable frame is 3.
    These 3 landmarks become constants in the problem, not variables to be optimized.

    The remaining landmarks are variables in the optimization.
    Let N be the total number of landmarks.
    The number of variable landmarks = Total Landmarks - Fixed Landmarks.

    The Schur complement strategy is used to marginalize these landmark *variables* to
    efficiently solve the optimization problem.

    Therefore, the maximum number of landmarks that can be marginalized is the number of
    landmarks that are treated as variables.
    """

    N_symbol = "N"
    fixed_landmarks = 3

    print("Explanation:")
    print("1. A bundle adjustment problem has a 7-degree-of-freedom gauge ambiguity (3 translation, 3 rotation, 1 scale).")
    print("2. To get a unique solution, this ambiguity must be removed by fixing a reference frame.")
    print("3. The minimum number of landmarks required to define a stable 3D reference frame is 3.")
    print("4. These 3 landmarks are treated as fixed constants, not as variables to be solved for.")
    print(f"5. The number of landmarks that remain as variables is: {N_symbol} - {fixed_landmarks}")
    print("6. The Schur complement is used to marginalize the landmark *variables*.")
    print("\nFinal Equation:")
    print(f"Maximum number of landmarks that can be marginalized = {N_symbol} - {fixed_landmarks}")

if __name__ == '__main__':
    solve_bundle_adjustment_question()