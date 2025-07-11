import numpy as np
from scipy.optimize import fsolve

def solve_repressilator_bifurcation():
    """
    This script finds the critical value of 'n' for the onset of oscillations
    in the given biochemical system by locating the Hopf bifurcation point.
    """
    # System parameters
    alpha = 100
    beta = 20

    print("Analyzing the system of equations:")
    print("dx/dt = alpha / (1 + z^n) - x - beta / (1 + x)")
    print("dy/dt = alpha / (1 + x^n) - y - beta / (1 + y)")
    print("dz/dt = alpha / (1 + y^n) - z - beta / (1 + z)")
    print(f"with alpha = {alpha} and beta = {beta}\n")

    def bifurcation_equations(variables):
        """
        Defines the system of two non-linear equations to find the critical
        point (x_star, n) for the Hopf bifurcation.
        
        Args:
            variables (list): A list containing [x_star, n].
        
        Returns:
            list: The residuals of the two equations ([eq1, eq2]).
        """
        x_star, n = variables
        
        # Guard against non-physical values (negative concentration or n)
        # that the solver might try during its search.
        if x_star <= 0 or n <= 0:
            return [1e6, 1e6]

        try:
            # Pre-calculate powers to avoid re-computation and handle potential overflow
            x_pow_n = np.power(x_star, n)
            
            # Check for overflow
            if np.isinf(x_pow_n):
                # For very large x**n, the term alpha/(1+x**n) -> 0
                eq1 = -x_star - beta / (1 + x_star)
                eq2 = -1 + beta / (1 + x_star)**2
            else:
                x_pow_n_minus_1 = np.power(x_star, n - 1)

                # Equation 1: The steady-state condition f(x*, n) = 0
                eq1 = alpha / (1 + x_pow_n) - x_star - beta / (1 + x_star)
                
                # Equation 2: The Hopf bifurcation condition, Re(lambda) = A - B/2 = 0
                # 'A' and 'B' are terms from the system's Jacobian matrix.
                A = -1 + beta / (1 + x_star)**2
                B = -alpha * n * x_pow_n_minus_1 / (1 + x_pow_n)**2
                eq2 = A - B / 2

        except (ValueError, OverflowError):
            # Return a large residual if a math error occurs
            return [1e6, 1e6]

        return [eq1, eq2]

    # Provide an initial guess for the solver: [x_star_guess, n_guess].
    # A good guess for x_star is around sqrt(beta).
    initial_guess = [np.sqrt(beta), 2.0]
    
    # Use fsolve to find the root [x_star_crit, n_crit] of the system
    solution, info, ier, msg = fsolve(bifurcation_equations, initial_guess, full_output=True)

    if ier == 1:
        x_star_crit, n_crit = solution
        print("A Hopf bifurcation occurs, leading to oscillations, when the parameter 'n' exceeds a critical value.")
        print("This critical value is found by numerically solving the steady-state and stability equations simultaneously.")
        print("\n--- Results ---")
        print(f"The steady-state concentration at the bifurcation point is x* = {x_star_crit:.4f}")
        print(f"The critical value for the parameter 'n' is n_crit = {n_crit:.4f}")
        
        # The final condition for oscillations, representing the final equation
        print("\nThe final equation describing the condition for the system to exhibit oscillations is:")
        print(f"n > {n_crit:.4f}")
        
        answer_string = f"n > {n_crit:.2f}"
        print(f"\n<<<{answer_string}>>>")
        
    else:
        print("The numerical solver failed to find a solution.")
        print("Solver message:", msg)
        print("\n<<<Solver failed>>>")

# Execute the main function
solve_repressilator_bifurcation()