import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

def solve_repressilator_problem():
    """
    Finds the critical value of n for oscillations in a repressilator system.
    """
    # Parameters from the problem description
    alpha = 100.0
    beta = 20.0

    # The system of ODEs describing the repressilator
    # dx/dt = alpha / (1 + z^n) - x - beta / (1 + x)
    # dy/dt = alpha / (1 + x^n) - y - beta / (1 + y)
    # dz/dt = alpha / (1 + y^n) - z - beta / (1 + z)
    def repressilator_odes(t, Y, n, alpha, beta):
        x, y, z = Y
        dxdt = alpha / (1 + z**n) - x - beta / (1 + x)
        dydt = alpha / (1 + x**n) - y - beta / (1 + y)
        dzdt = alpha / (1 + y**n) - z - beta / (1 + z)
        return [dxdt, dydt, dzdt]

    # Helper function to check if the system oscillates for a given n
    def is_oscillating(n, alpha, beta):
        """
        Simulates the system for a given n and returns True if it oscillates.
        """
        # 1. Find the symmetric steady state (x=y=z=x_ss) by finding the root of:
        # alpha / (1 + x_ss^n) - x_ss - beta / (1 + x_ss) = 0
        def steady_state_eq(x):
            # Ensure x is positive for physical concentrations
            if x <= 0:
                return float('inf')
            return alpha / (1 + x**n) - x - beta / (1 + x)

        try:
            # Use fsolve to find the root, starting with a reasonable guess
            x_ss, = fsolve(steady_state_eq, x0=2.0)
            if x_ss <= 0:
                return False
        except:
            # If root finding fails, assume no stable oscillation for this n
            return False

        # 2. Simulate the ODE system from a perturbed initial condition
        y0 = [x_ss * 1.05, x_ss, x_ss]  # 5% perturbation on x
        t_span = [0, 500]
        
        # Use solve_ivp for numerical integration
        sol = solve_ivp(
            fun=repressilator_odes,
            t_span=t_span,
            y0=y0,
            method='RK45',
            args=(n, alpha, beta),
            dense_output=True
        )

        if not sol.success:
            return False

        # 3. Analyze the solution to detect oscillations
        # Evaluate solution on a fine grid to analyze the tail end
        t_eval = np.linspace(t_span[0], t_span[1], 1000)
        sol_y = sol.sol(t_eval)

        # Check the last 25% of the simulation for oscillations
        transient_index = int(0.75 * len(t_eval))
        x_series_tail = sol_y[0, transient_index:]
        
        # Calculate amplitude and mean of the tail
        amplitude = np.max(x_series_tail) - np.min(x_series_tail)
        mean_val = np.mean(x_series_tail)

        # Consider it oscillating if relative amplitude is > 1%
        if mean_val > 1e-6 and (amplitude / mean_val) > 0.01:
            return True
        else:
            return False

    # Use binary search to find the critical value of n
    low_n = 1.0
    high_n = 5.0
    precision = 0.01

    # Check that oscillations exist within the search range
    if not is_oscillating(high_n, alpha, beta):
        print(f"Warning: No oscillations detected for n up to {high_n}.")
        print("The critical value might be higher, or the system does not oscillate.")
        return

    while (high_n - low_n) > precision:
        mid_n = (low_n + high_n) / 2.0
        if is_oscillating(mid_n, alpha, beta):
            # If it oscillates, the critical value is lower or equal
            high_n = mid_n
        else:
            # If stable, the critical value is higher
            low_n = mid_n
    
    critical_n = high_n

    # Print the final result including the numbers from the "final equation"
    print("For the system of equations with alpha = 100 and beta = 20:")
    print("dx/dt = 100 / (1 + z^n) - x - 20 / (1 + x)")
    print("dy/dt = 100 / (1 + x^n) - y - 20 / (1 + y)")
    print("dz/dt = 100 / (1 + y^n) - z - 20 / (1 + z)")
    print("\n" + "="*50)
    print(f"The system exhibits oscillations for values of n > {critical_n:.2f}")
    print("="*50)

if __name__ == '__main__':
    solve_repressilator_problem()