import sys

def solve_ba_marginalization(N_symbolic="N", M_symbolic="M"):
    """
    Explains the reasoning to find the maximum number of landmarks that can be
    marginalized in a bundle adjustment problem.

    Args:
        N_symbolic (str): The symbol for the total number of landmarks.
        M_symbolic (str): The symbol for the total number of cameras.
    """
    print(f"Problem: Given {N_symbolic} landmarks and {M_symbolic} cameras in bundle adjustment,")
    print("what is the maximum number of landmarks that can be marginalized using the Schur complement,")
    print("while ensuring the problem remains solvable?")
    print("-" * 70)

    print("Step 1: Understand Gauge Freedom in Bundle Adjustment")
    print("The bundle adjustment solution is ambiguous up to a similarity transformation (3D rotation, 3D translation, and scale).")
    num_gauge_dof = 7
    print(f"This results in a {num_gauge_dof}-degree-of-freedom (DOF) gauge freedom that must be fixed for a unique solution.")
    print("\nStep 2: Interpret 'Solvable'")
    print("A 'solvable' problem in this context means one where a unique metric reconstruction can be found.")
    print("This requires fixing the gauge freedom.")

    print("\nStep 3: Define a Stable Reference Frame")
    print("A stable 3D reference frame can be defined by fixing the positions of a set of non-collinear landmarks.")
    num_landmarks_for_frame = 3
    print(f"A minimum of {num_landmarks_for_frame} landmarks are needed to unambiguously fix the 6-DOF for rotation and translation.")
    
    print("\nStep 4: Fix the Scale")
    print("The 7th DOF (scale) is fixed by setting a unit distance within the reference frame, using the same landmarks.")
    print(f"Therefore, a total of {num_landmarks_for_frame} landmarks are used to fix all {num_gauge_dof} DOFs.")

    print("\nStep 5: Relate to the Schur Complement and Marginalization")
    print("Marginalization eliminates variables from the optimization system.")
    print("The landmarks used to define the fixed reference frame are no longer free variables.")
    print("As they are not part of the optimization variables, they cannot be marginalized.")

    print("\nStep 6: Form the Final Equation")
    print("The maximum number of landmarks that can be marginalized is the total number of landmarks")
    print("minus the minimum number of landmarks required to establish the fixed reference frame.")
    
    print("\nFinal Equation:")
    # Printing each 'number' in the final equation as requested
    print(f"Max Marginalizable Landmarks = {N_symbolic} - {num_landmarks_for_frame}")

# Execute the explanation
solve_ba_marginalization()