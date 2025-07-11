def solve():
    """
    Calculates the maximum number of landmarks that can be marginalized in a
    bundle adjustment problem.
    """
    # In bundle adjustment, the system of cameras and landmarks has a 7-DOF
    # ambiguity (3 for translation, 3 for rotation, 1 for scale).
    # To make the problem solvable, this ambiguity must be removed by fixing
    # a coordinate system, which is known as fixing the gauge.

    # A fundamental way to fix the gauge in 3D space is to define a reference
    # frame using a minimal number of landmarks.
    # A minimum of 3 non-collinear landmarks are required for this.
    # - Landmark 1 fixes the origin (3 DOFs).
    # - Landmark 2 fixes an axis and scale (3 DOFs).
    # - Landmark 3 fixes the final rotation (1 DOF).
    num_fixed_landmarks = 3

    # The Schur complement strategy works by marginalizing (eliminating) landmark
    # variables to create a smaller system for the camera variables.
    # The landmarks used to fix the gauge are no longer variables to be solved for;
    # they are constants. Therefore, they cannot be marginalized.

    # If we have N total landmarks, the maximum number that can be treated as
    # variables and therefore marginalized is the total number minus the number
    # required to fix the gauge.
    
    # We use a symbolic representation for N in the explanation.
    # For the final printout, we'll use an example value.
    example_n = 20

    max_marginalized = example_n - num_fixed_landmarks

    print("The bundle adjustment problem has a 7-DOF gauge ambiguity that must be fixed for a unique solution.")
    print(f"To fix this gauge, a minimum of {num_fixed_landmarks} landmarks are required to define a fixed coordinate system.")
    print("These fixed landmarks are no longer variables and cannot be marginalized.")
    print("Given N total landmarks, the maximum number that can be marginalized is N - 3.")
    print("\nFor example, if N is 20:")
    print(f"Equation: Total Landmarks - Fixed Landmarks = Max Marginalized Landmarks")
    print(f"Numbers:  {example_n} - {num_fixed_landmarks} = {max_marginalized}")

solve()