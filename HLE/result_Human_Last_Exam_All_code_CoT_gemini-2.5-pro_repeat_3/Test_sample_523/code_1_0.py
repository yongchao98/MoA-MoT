import numpy as np
from scipy.optimize import fsolve, brentq

def find_oscillation_boundary():
    """
    This script finds the critical value of the parameter 'n' for which the
    biochemical system starts to exhibit oscillations.
    """
    # Define system parameters
    alpha = 100.0
    beta = 20.0

    print(f"System parameters: alpha = {alpha}, beta = {beta}\n")

    def bifurcation_condition(n):
        """
        Calculates the value of the Hopf bifurcation condition, Re(lambda) = A - B/2.
        The system oscillates when this value is > 0.
        The bifurcation occurs when this value is = 0.
        This function is designed to be used with a root finder to find the
        critical value of n.
        """
        # Define the steady-state equation to find 's' for the given 'n'
        def steady_state_eq(s_vec):
            s = s_vec[0]
            # s must be a positive concentration
            if s <= 0:
                return 1e9 # Return a large error if s is not positive
            return alpha / (1 + s**n) - s - beta / (1 + s)

        # Numerically solve for the steady state 's'
        # An initial guess of s=2.0 is robust for this problem.
        s_star_array, info, ier, msg = fsolve(steady_state_eq, x0=[2.0], full_output=True)
        if ier != 1:
            # fsolve failed to converge, return a large value to indicate this 'n' is not the root
            return 1e9
        s_star = s_star_array[0]

        # Calculate the terms A and B of the Jacobian eigenvalues
        A = -1.0 + beta / (1.0 + s_star)**2
        B = -alpha * n * s_star**(n - 1.0) / (1.0 + s_star**n)**2

        # The bifurcation condition
        return A - B / 2.0

    try:
        # Find the critical value of n by finding the root of the bifurcation_condition function.
        # We search in the interval [1, 10] because we know the function changes sign in this range.
        n_critical = brentq(bifurcation_condition, a=1.0, b=10.0)

        print(f"The system exhibits oscillations for n > n_critical.")
        print(f"The critical value n_critical is found by numerically solving the bifurcation condition A - B/2 = 0.")

        # For the final output, let's re-calculate the values at the critical point
        # to display the final equation with numbers.
        def ss_eq_final(s_vec):
            s = s_vec[0]
            return alpha / (1 + s**n_critical) - s - beta / (1 + s)
        s_final_array = fsolve(ss_eq_final, x0=[2.0])
        s_final = s_final_array[0]

        A_final = -1.0 + beta / (1.0 + s_final)**2
        B_final = -alpha * n_critical * s_final**(n_critical - 1.0) / (1.0 + s_final**n_critical)**2
        
        print(f"\n--- Calculation at the Bifurcation Point ---")
        print(f"Critical value of n found: n_critical = {n_critical:.4f}")
        print(f"Steady state concentration at this point: s = {s_final:.4f}")
        
        print("\nThe final equation for the bifurcation is A - B/2 = 0, where:")
        print(f"A = -1 + {beta:.1f} / (1 + {s_final:.4f})^2 = {A_final:.4f}")
        print(f"B = -{alpha:.1f} * {n_critical:.4f} * {s_final:.4f}^({n_critical:.4f}-1) / (1 + {s_final:.4f}^{n_critical:.4f})^2 = {B_final:.4f}")
        
        print("\nSubstituting these values into the bifurcation equation:")
        print(f"A - B/2 = {A_final:.4f} - ({B_final:.4f}) / 2 = {A_final - B_final / 2.0:.6f}")
        print("(The result is approximately zero, as expected at the bifurcation point.)")

        print(f"\nCONCLUSION: The system exhibits oscillations for values of n > {n_critical:.4f}.")
        
        # Storing the final answer for the specified format
        global answer
        answer = n_critical

    except ValueError:
        print("Could not find the critical value of n in the specified range.")
        answer = "Error"

find_oscillation_boundary()
print(f"\n<<<{answer}>>>")