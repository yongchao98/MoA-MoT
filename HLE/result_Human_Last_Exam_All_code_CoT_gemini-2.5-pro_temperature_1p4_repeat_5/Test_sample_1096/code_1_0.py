import math

def solve_optimal_steps():
    """
    Calculates the optimal two-step gradient descent learning rates (gamma_1, gamma_2)
    for a given condition number kappa.
    """
    # We consider an M-smooth and mu-strongly convex function.
    # The problem assumes mu=1 and M=kappa. Let's use kappa=10 for a numerical example.
    kappa = 10.0

    print(f"The problem is to find the best choice of (gamma_1, gamma_2) for a 2-step gradient descent.")
    print(f"This choice minimizes the convergence rate uniformly over the class of functions.")
    print(f"The analysis reduces to a minimax polynomial problem on the interval of eigenvalues [1, kappa].")
    print(f"Assuming a condition number kappa = {kappa}:\n")
    
    # The optimal step sizes are derived from the solution to the minimax problem.
    # The derived formulas in terms of kappa are:
    # gamma_1 = (4*(kappa+1) - 2*sqrt(2)*(kappa-1)) / (kappa^2 + 6*kappa + 1)
    # gamma_2 = (4*(kappa+1) + 2*sqrt(2)*(kappa-1)) / (kappa^2 + 6*kappa + 1)
    # The order of gamma_1 and gamma_2 can be swapped.
    
    # Numerator components
    term1_num = 4 * (kappa + 1)
    term2_num = 2 * math.sqrt(2) * (kappa - 1)
    
    # Denominator
    denominator = kappa**2 + 6 * kappa + 1
    
    # Calculate gamma_1 and gamma_2
    gamma1 = (term1_num - term2_num) / denominator
    gamma2 = (term1_num + term2_num) / denominator
    
    print("The symbolic formulas for the optimal step sizes are:")
    print("gamma_1 = (4*(kappa+1) - 2*sqrt(2)*(kappa-1)) / (kappa^2 + 6*kappa + 1)")
    print("gamma_2 = (4*(kappa+1) + 2*sqrt(2)*(kappa-1)) / (kappa^2 + 6*kappa + 1)\n")
    
    print(f"For kappa = {kappa}, the numerical values are:")
    print(f"gamma_1 = {gamma1}")
    print(f"gamma_2 = {gamma2}\n")
    
    print("--- Breakdown of the Final Equation ---")
    print(f"Value of kappa: {kappa}")
    print(f"Value of mu: 1.0")
    print(f"Numerator Term 1 (4*(kappa+1)): {term1_num}")
    print(f"Numerator Term 2 (2*sqrt(2)*(kappa-1)): {term2_num}")
    print(f"Denominator (kappa^2 + 6*kappa + 1): {denominator}")
    
    print(f"\nFinal equation for gamma_1:")
    print(f"{gamma1} = ({term1_num:.4f} - {term2_num:.4f}) / {denominator:.4f}")
    
    print(f"\nFinal equation for gamma_2:")
    print(f"{gamma2} = ({term1_num:.4f} + {term2_num:.4f}) / {denominator:.4f}")

if __name__ == '__main__':
    solve_optimal_steps()