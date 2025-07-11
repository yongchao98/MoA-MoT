import numpy as np
from scipy.optimize import fsolve

def find_oscillation_threshold():
    """
    This script calculates the critical value of the parameter 'n' for which the
    given biochemical system starts to exhibit oscillations.
    """
    # System parameters
    alpha = 100.0
    beta = 20.0

    # --- Step 1: Define functions to find the steady state ---
    def steady_state_equation(x, n, alpha, beta):
        """Equation for the steady state x_ss. The root of this is x_ss."""
        return alpha / (1 + x**n) - x - beta / (1 + x)

    def get_steady_state(n, alpha, beta):
        """Finds the steady state x_ss for a given n using a numerical solver."""
        # Initial guess for x_ss. Based on analysis, the root is usually > 1.
        initial_guess_x = 5.0
        x_ss, = fsolve(steady_state_equation, initial_guess_x, args=(n, alpha, beta))
        return x_ss

    # --- Step 2: Define the bifurcation condition ---
    def bifurcation_condition(n, alpha, beta):
        """
        This function represents the real part of the complex eigenvalues (A + B/2).
        A Hopf bifurcation occurs when this function is zero.
        """
        if n <= 0:
            return np.inf  # n must be positive
        
        # Get the steady state for the given n
        x_ss = get_steady_state(n, alpha, beta)
        
        if x_ss <= 0:
            return np.inf # Steady state must be physically meaningful (positive)

        # Calculate Jacobian elements A and B
        A = -1.0 + beta / (1.0 + x_ss)**2
        B = alpha * n * x_ss**(n-1) / (1.0 + x_ss**n)**2
        
        # Return the value of the stability condition
        return A + B / 2.0

    # --- Step 3: Solve for the critical value of n ---
    # Based on manual checks, the root for n is between 1 and 2.
    initial_guess_n = 1.5
    n_critical, = fsolve(bifurcation_condition, initial_guess_n, args=(alpha, beta))

    # --- Step 4: Output the results and verification ---
    print(f"For the system with parameters alpha = {alpha} and beta = {beta}:")
    print("Oscillations occur when the real part of the complex eigenvalues becomes positive.")
    print("This happens at a critical value of 'n' where the condition 'A + B/2 = 0' is met.")
    
    # Calculate values at the critical point for display
    x_ss_critical = get_steady_state(n_critical, alpha, beta)
    A_critical = -1.0 + beta / (1.0 + x_ss_critical)**2
    B_critical = alpha * n_critical * x_ss_critical**(n_critical-1) / (1.0 + x_ss_critical**n_critical)**2
    
    print("\n--- Results at the Bifurcation Point ---")
    print(f"The critical value for n is: {n_critical:.4f}")
    print(f"The steady state concentration (x_ss) at this point is: {x_ss_critical:.4f}")
    
    print("\nVerification of the final equation 'A + B/2 = 0':")
    # Print the equation with all numbers substituted
    print(f"A = -1 + {beta} / (1 + {x_ss_critical:.4f})^2 = {A_critical:.4f}")
    print(f"B = {alpha} * {n_critical:.4f} * {x_ss_critical:.4f}^({n_critical-1:.4f}) / (1 + {x_ss_critical:.4f}^{n_critical:.4f})^2 = {B_critical:.4f}")
    
    result = A_critical + B_critical / 2.0
    print(f"A + B/2 = {A_critical:.4f} + {B_critical:.4f} / 2 = {result:.4f}")
    
    print("\nConclusion:")
    print(f"The system exhibits oscillations for n > {n_critical:.4f}.")
    print("The final answer is the critical value where oscillations begin.")
    
    return n_critical

# Execute the function and print the final answer in the required format
critical_n_value = find_oscillation_threshold()
print(f"\n<<<>>>\n{critical_n_value}\n<<<>>>")
