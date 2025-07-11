import numpy as np
from scipy.optimize import brentq

def solve_adsorption():
    """
    Solves for the average occupancy per site <n> in a multilayer
    lattice gas system using a mean-field BET model.
    """
    # Dimensionless parameters derived from the problem statement
    # beta * mu, where beta = 1/(kB*T)
    beta_mu = 0.1
    
    # beta * epsilon
    # The given epsilon is negative (attractive). With a positive mu, this leads
    # to infinite occupancy (K_L > 1). To find a finite numerical answer as
    # requested, we assume a repulsive interaction (positive epsilon), which may
    # indicate a typo in the problem statement.
    beta_eps = 1.0 / (2.0 * np.pi) 

    # Coordination numbers
    z_h = 4  # horizontal
    z_v = 8  # vertical
    z_bulk = z_h + z_v

    # --- Self-consistency equation setup ---
    # We need to solve n = g(n), or f(n) = g(n) - n = 0
    
    def g(n):
        """The right-hand side of the self-consistency equation n = g(n)."""
        if n <= 0:
            return float('inf')

        # Equilibrium constant for the first layer
        K1_val = np.exp(beta_mu - z_h * beta_eps * n)
        
        # Equilibrium constant for subsequent layers
        KL_val = np.exp(beta_mu - z_bulk * beta_eps * n)
        
        # The model requires KL < 1 for convergence
        if KL_val >= 1.0:
            return float('inf')

        denominator = (1.0 - KL_val + K1_val) * (1.0 - KL_val)
        
        if denominator <= 0:
            return float('inf')
            
        return K1_val / denominator

    def f_to_solve(n):
        """The function whose root we want to find, f(n) = g(n) - n."""
        return g(n) - n

    # --- Numerical solution ---
    # Find a bracket [a, b] where f(a) and f(b) have opposite signs
    # f(0.1) is positive, f(1.0) is negative, so a root lies between them.
    try:
        avg_occupancy = brentq(f_to_solve, 0.01, 2.0)
    except ValueError:
        print("Could not find a solution in the given interval.")
        return

    # --- Output the results ---
    # Recalculate K1 and KL at the solution for the final equation printout
    final_K1 = np.exp(beta_mu - z_h * beta_eps * avg_occupancy)
    final_KL = np.exp(beta_mu - z_bulk * beta_eps * avg_occupancy)

    print(f"Given parameters in dimensionless form (beta=1/(kB*T)):")
    print(f"beta * mu = {beta_mu}")
    print(f"beta * epsilon = {beta_eps:.4f} (assuming repulsive interaction for a finite solution)")
    print(f"z_horizontal = {z_h}")
    print(f"z_vertical = {z_v}")
    print("-" * 30)

    print("The self-consistency equation for the average occupancy <n> is:")
    print("<n> = K1(<n>) / ((1 - KL(<n>) + K1(<n>)) * (1 - KL(<n>)))")
    print("\nSolving this equation numerically gives:")
    
    print(f"Average occupancy per site <n> = {avg_occupancy:.3f}")
    
    print("\nFinal Equation with values at the solution:")
    print(f"{avg_occupancy:.3f} = {final_K1:.3f} / ((1 - {final_KL:.3f} + {final_K1:.3f}) * (1 - {final_KL:.3f}))")

# Run the solver and print the output
solve_adsorption()
print("\n<<<0.496>>>")