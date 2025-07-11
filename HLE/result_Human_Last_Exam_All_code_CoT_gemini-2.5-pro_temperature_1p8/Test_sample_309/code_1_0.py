def calculate_max_marginalized_landmarks():
    """
    This function explains the reasoning and derives the formula for the maximum
    number of landmarks that can be marginalized in a bundle adjustment problem
    while ensuring solvability.
    """

    print("Step-by-step reasoning to determine the maximum number of marginalizable landmarks:")
    print("-" * 70)

    print("\n[Step 1] The Problem: We have N landmarks and M cameras in a Bundle Adjustment (BA) setup.")
    print("The goal is to use the Schur Complement to make the optimization efficient by marginalizing landmarks.")

    print("\n[Step 2] The Core Constraint: The problem must remain 'solvable'.")
    print("In 3D reconstruction, 'solvable' means the resulting structure is rigid and not degenerate.")
    print("The solution should be unique up to a global translation, rotation, and scale (7 degrees of freedom).")

    print("\n[Step 3] The Requirement for a Rigid Structure: A stable geometric skeleton.")
    print("To define a rigid 3D structure, we need a minimum set of anchor points:")
    print(" - 1 point is not enough (no shape).")
    print(" - 2 points define a line, but the rest of the scene can rotate around this line (not rigid in 3D).")
    print(" - 3 non-collinear points define a rigid triangle. This is the minimal skeleton for a stable 3D scene.")

    print("\n[Step 4] Connecting Rigidity to Marginalization:")
    print("The landmarks that form this essential 3-point skeleton must be explicitly kept in the optimization state to act as anchors.")
    print("Therefore, these anchor landmarks CANNOT be marginalized.")
    
    min_anchor_landmarks = 3
    
    print(f"\n[Step 5] The Calculation:")
    print(f"To ensure a solvable system, a minimum of {min_anchor_landmarks} landmarks must be kept.")
    print("If we have N total landmarks, the maximum number we can afford to marginalize is the total minus this required minimum.")

    print("\nFinal Equation:")
    print("Max Marginalized Landmarks = Total Landmarks (N) - Minimum Anchor Landmarks")
    print(f"Resulting Formula = N - {min_anchor_landmarks}")

calculate_max_marginalized_landmarks()