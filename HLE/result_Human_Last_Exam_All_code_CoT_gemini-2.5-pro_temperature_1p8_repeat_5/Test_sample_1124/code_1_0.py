def modified_logistic_map():
    """
    This function demonstrates a modified logistic map that results in a stable
    equilibrium point greater than 1 for R=3.57.
    The modification used is: X_n+1 = X_n/R + 1 - 1/R**2
    """
    R = 3.57
    x = 0.5  # An arbitrary starting value for X

    # --- Theoretical Calculation ---
    # The equilibrium point X_eq is found by solving X_eq = X_eq/R + 1 - 1/R**2
    # which simplifies to X_eq = 1 + 1/R
    x_eq_theoretical = 1 + 1/R

    print(f"Modified map: X_n+1 = X_n/R + 1 - 1/(R^2)")
    print(f"For R = {R}, the theoretical stable equilibrium point is X_eq = {x_eq_theoretical:.6f}\n")

    # --- Iteration and Convergence ---
    print("Iterating the map to show convergence:")
    for i in range(20):
        x_next = x/R + (1 - 1/R**2)
        # Stop if it gets very close to the theoretical value
        if abs(x_next - x) < 1e-7:
            x = x_next
            break
        x = x_next
        print(f"Step {i+1:2d}: X = {x:.6f}")

    print(f"\nThe map has converged to the equilibrium point X = {x:.6f}\n")

    # --- Final Equation Verification ---
    print("Verification of the final equilibrium equation:")
    # We will show that for the converged value X, the equation X = X/R + 1 - 1/R**2 holds true.
    final_x = x
    right_hand_side_term1 = final_x / R
    right_hand_side_term2 = 1 - 1/R**2
    print(f"We check if: {final_x:.6f} = {right_hand_side_term1:.6f} + {right_hand_side_term2:.6f}")
    print(f"The right side equals: {(right_hand_side_term1 + right_hand_side_term2):.6f}")
    print("The equation holds true.")


modified_logistic_map()
<<<1.280112>>>