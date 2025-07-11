import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import warnings

# Suppress warnings that can occur during root finding for clarity
warnings.filterwarnings('ignore', 'The iteration is not making good progress')

def find_oscillation_threshold():
    """
    This script finds the critical value of the parameter 'n' for which the
    given dynamical system starts to exhibit oscillations.
    """

    # System parameters from the problem description
    alpha = 100.0
    beta = 20.0

    def system_ode(t, y, n):
        """Defines the system of three coupled ODEs."""
        x, y, z = y
        dxdt = alpha / (1 + z**n) - x - beta / (1 + x)
        dydt = alpha / (1 + x**n) - y - beta / (1 + y)
        dzdt = alpha / (1 + y**n) - z - beta / (1 + z)
        return [dxdt, dydt, dzdt]

    def find_fixed_point(n):
        """Numerically finds the symmetric fixed point u for a given n."""
        # The fixed point equation is: u + beta/(1+u) - alpha/(1+u^n) = 0
        def fp_eq(u):
            # The physical concentrations must be non-negative
            if u[0] < 0:
                return 1e6
            return u[0] + beta / (1 + u[0]) - alpha / (1 + u[0]**n)

        # fsolve finds the root of the equation. We provide an initial guess (e.g., 5.0)
        u_star, _, exit_flag, _ = fsolve(fp_eq, x0=[5.0], full_output=True)
        if exit_flag == 1:
            return u_star[0]
        else:
            return None # Indicate that a stable fixed point was not found

    def check_for_oscillations(n):
        """
        For a given n, simulates the system and returns True if it oscillates,
        False otherwise.
        """
        u_star = find_fixed_point(n)
        if u_star is None:
            # If a simple fixed point isn't found, behavior might be complex.
            # Assume it could be oscillatory for the purpose of this search.
            return True

        # Initial condition: Perturb the system from its fixed point
        y0 = [u_star * 1.05, u_star, u_star]

        # Time span for the simulation
        t_span = [0, 500]
        
        # Solve the ODE
        sol = solve_ivp(
            system_ode, 
            t_span, 
            y0, 
            args=(n,), 
            method='RK45', 
            dense_output=True
        )

        # Analyze the tail end of the simulation for oscillations
        if sol.status == 0:
            t_eval = np.linspace(t_span[1] * 0.75, t_span[1], 500)
            y_eval = sol.sol(t_eval)
            x_late = y_eval[0]
            
            # Oscillation detection: check if the amplitude is significant
            amplitude = np.max(x_late) - np.min(x_late)
            # A small threshold distinguishes oscillations from convergence to a fixed point
            if amplitude > 1.0:
                return True
        
        return False

    # Search for the critical value of n by testing values from 1.0 upwards
    print("Searching for the critical value of n...")
    critical_n = -1
    for n_test in np.arange(1.0, 2.5, 0.01):
        if check_for_oscillations(n_test):
            critical_n = n_test
            break
            
    if critical_n > 0:
        # The final equation describes the condition on n for oscillations
        print("\n--- Final Result ---")
        print("The dynamical system is defined by:")
        print("dx/dt = alpha / (1 + z^n) - x - beta / (1 + x)")
        print("dy/dt = alpha / (1 + x^n) - y - beta / (1 + y)")
        print("dz/dt = alpha / (1 + y^n) - z - beta / (1 + z)")
        print("\nWith the given parameter values:")
        print(f"alpha = {int(alpha)}")
        print(f"beta = {int(beta)}")
        print("\nThe system exhibits oscillations for values of n where:")
        # We output the final equation/condition as requested
        print(f"n > {critical_n - 0.01:.2f}")
    else:
        print("Could not determine the threshold for oscillations in the tested range.")
        
    return critical_n

if __name__ == '__main__':
    find_oscillation_threshold()