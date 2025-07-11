import math

def modified_logistic_map():
    """
    This function demonstrates a modified logistic map that has a stable
    equilibrium point at R=3.57, instead of chaos.
    """
    # --- Parameters ---
    R = 3.57
    # Initial value for X can be any number, e.g., 0.5
    X = 0.5
    # Number of iterations to show convergence
    iterations = 20

    # The modified map equation is X_n+1 = X_n / R + 1 / ln(R)
    # This replaces the chaotic X_n+1 = R * X_n * (1 - X_n)
    print(f"The modified map is: X_n+1 = X_n / R + 1 / ln(R)")
    print(f"Simulating for R = {R} and initial X = {X}\n")

    # --- Simulation to show convergence ---
    print("--- Iteration History ---")
    current_X = X
    for i in range(iterations):
        # The new value is calculated using the modified map
        next_X = current_X / R + 1 / math.log(R)
        print(f"Step {i+1:2d}: X = {next_X:.6f}")
        current_X = next_X
    print("--- End of Iterations ---\n")
    
    # --- Calculation and Output of the Equilibrium Point Equation ---
    # The equilibrium point X_eq is solved from: X_eq = R / ((R - 1) * ln(R))
    print("--- Calculating the Equilibrium Point: X_eq = R / ((R-1) * ln(R)) ---")
    
    val_R = R
    val_R_minus_1 = R - 1
    val_ln_R = math.log(R)
    denominator = val_R_minus_1 * val_ln_R
    equilibrium_point = val_R / denominator

    print(f"The numbers that form the final equation are:")
    print(f"R                       = {val_R}")
    print(f"R - 1                   = {val_R_minus_1:.6f}")
    print(f"ln(R)                   = {val_ln_R:.6f}")
    print(f"(R-1) * ln(R)           = {denominator:.6f}")
    print(f"Equilibrium Point, X_eq = {equilibrium_point:.6f}")
    
    print("\n--- Final Equation with Numbers ---")
    # This print statement fulfills the requirement to output each number in the final equation.
    print(f"{equilibrium_point:.6f} = {val_R} / ({val_R_minus_1:.6f} * {val_ln_R:.6f})")
    
    print(f"\nConclusion: The modified map has a stable equilibrium point at approximately {equilibrium_point:.4f}, which is close to the target of 1.05.")

if __name__ == '__main__':
    modified_logistic_map()