import numpy as np
from scipy.optimize import fsolve

def find_oscillatory_regime():
    """
    This function calculates the critical value of the parameter 'n' for which the
    given biochemical system starts to exhibit oscillations.
    """
    # System parameters as given in the problem
    alpha = 100.0
    beta = 20.0

    # We define a function that represents the bifurcation condition.
    # The root of this function will be the critical value of 'n'.
    # This function internally solves for the steady state for any given 'n'.
    def bifurcation_condition(n_array):
        """
        Calculates the value of the stability condition (Re(lambda_complex))
        for a given n. The root of this function corresponds to the Hopf bifurcation.
        fsolve requires the function to take a 1-D array.
        """
        n = n_array[0]

        # Step 1: Find the steady state (x_ss) for the given n.
        # We need to solve the equation: x + beta/(1+x) - alpha/(1+x^n) = 0
        def steady_state_equation(x, n_val):
            if x <= 0:  # Concentration must be positive
                return np.inf
            return x + beta / (1 + x) - alpha / (1 + x**n_val)

        # Use a numerical solver (fsolve) to find x_ss.
        # An initial guess of 3.0 is found to be robust for the relevant range of n.
        x_ss_solution = fsolve(steady_state_equation, x0=3.0, args=(n,))
        x_ss = x_ss_solution[0]

        # Step 2: Calculate the terms 'A' and 'B' from the Jacobian matrix at the steady state.
        # A = J_ii = -1 + beta / (1 + x_ss)^2
        # B = J_ij (for coupled terms) = -alpha * n * x_ss^(n-1) / (1 + x_ss^n)^2
        A = -1.0 + beta / (1.0 + x_ss)**2

        # To avoid numerical overflow for large n and x > 1, we use a stable formulation for B.
        if x_ss > 1:
            B = -alpha * n * (x_ss**(-n - 1)) / (x_ss**(-n) + 1)**2
        else:
            B = -alpha * n * x_ss**(n - 1) / (1 + x_ss**n)**2
        
        # Step 3: The Hopf bifurcation occurs when the real part of the complex eigenvalues is zero.
        # The eigenvalues are lambda_1 = A+B and lambda_2,3 = A - B/2 ± i*sqrt(3)/2*B.
        # The condition for bifurcation is Re(lambda_2,3) = A - B/2 = 0.
        return A - B / 2.0

    # Step 4: Find the root of the bifurcation_condition function to get the critical n.
    # From preliminary analysis, the root is known to be between 2 and 3.
    # We use an initial guess of 2.8 for the solver.
    initial_guess_n = 2.8
    critical_n_solution = fsolve(bifurcation_condition, x0=initial_guess_n)
    critical_n = critical_n_solution[0]

    # Print the final result in a clear, descriptive format.
    print("The system exhibits oscillations when the parameter n is greater than a critical value.")
    print(f"This critical value, n_c, is found by solving the bifurcation condition numerically.")
    print("\nThe condition for oscillation is:")
    print(f"n > {critical_n:.4f}")

if __name__ == "__main__":
    find_oscillatory_regime()