import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar

def solve_critical_potential():
    """
    Calculates the critical chemical potential for scalar condensation in a holographic
    D3/D7 model with an Einstein-Gauss-Bonnet background.
    """
    # --- Model Parameters ---
    # Gauss-Bonnet coupling
    LAMBDA_GB = 0.1
    # Mass squared of the scalar field dual to the <bar{q}q> operator (dimension Delta=3)
    M_SQUARED = -3

    # --- EGB Metric Functions ---
    # Metric function f(z) where z = r_h / r. We set the horizon radius r_h = 1.
    def f(z, lgb):
        # Add a small epsilon to z to avoid division by zero at z=0
        z = z + 1e-12
        if lgb == 0:
            return 1 - z**4
        # Handle potential sqrt of negative number due to precision errors
        inner_sqrt_arg = 1 - 4 * lgb * (1 - z**4)
        if inner_sqrt_arg < 0:
            inner_sqrt_arg = 0
        return (1 / (2 * lgb * z**2)) * (1 - np.sqrt(inner_sqrt_arg))

    # Derivative of f(z) with respect to z
    def f_prime(z, lgb):
        z_safe = z + 1e-12
        if lgb == 0:
            return -4 * z_safe**3
        
        inner_sqrt_arg = 1 - 4 * lgb * (1 - z_safe**4)
        # Handle numerical stability near horizon z=1
        if abs(z-1.0) < 1e-7:
            return -4.0
        if inner_sqrt_arg <= 0:
            inner_sqrt_arg = 1e-12

        S = np.sqrt(inner_sqrt_arg)
        S_prime = 8 * lgb * z_safe**3 / S
        
        term1 = -(1 - S) / (lgb * z_safe**3)
        term2 = - (1 / (2 * lgb * z_safe**2)) * S_prime
        return term1 + term2

    # --- The Differential Equation System ---
    # We solve for Y = [psi, d(psi)/dz]
    def ode_system(z, Y, mu):
        psi, dpsi = Y
        
        f_val = f(z, LAMBDA_GB)
        f_prime_val = f_prime(z, LAMBDA_GB)

        # Avoid division by zero
        if abs(f_val) < 1e-9: f_val = 1e-9 if f_val > 0 else -1e-9
        if abs(z) < 1e-9: z = 1e-9
        
        # This ODE form ensures correct boundary behavior for a scalar with m^2=-3
        # psi'' + [f'/f - 3/z] psi' + [mu^2(1-z^2)^2/f^2 - m^2/(z^2 f)] psi = 0
        P_coeff = f_prime_val / f_val - 3 / z
        Q_coeff = (mu**2 * (1 - z**2)**2 / f_val**2) - (M_SQUARED / (z**2 * f_val))
        
        ddpsi = -P_coeff * dpsi - Q_coeff * psi
        return [dpsi, ddpsi]

    # --- Shooting Method Implementation ---
    # Objective function for the root finder. We want to find mu such that the
    # "source" term B in psi(z) = A*z + B*z^3 vanishes at the boundary z=0.
    def objective_function(mu):
        z_start = 1.0 - 1e-6 # Start just outside the horizon
        z_end = 1e-5       # End close to the boundary

        # Set initial conditions based on regularity at the horizon
        # For this ODE, psi'(1) = (m^2/4)*psi(1) = -3/4*psi(1). We set psi(1)=1.
        psi_0 = 1.0 - 0.75 * (1.0 - z_start)
        dpsi_0 = -0.75
        
        sol = solve_ivp(
            ode_system, [z_start, z_end], [psi_0, dpsi_0], args=(mu,),
            dense_output=True, method='RK45', atol=1e-8, rtol=1e-8
        )
        
        z_final = sol.t[-1]
        psi_final, dpsi_final = sol.y[:, -1]
        
        # From psi(z) = A*z + B*z^3 and psi'(z) = A + 3*B*z^2, we solve for B
        # We seek the value of mu that makes B = 0.
        B = (z_final * dpsi_final - psi_final) / (2 * z_final**3)
        return B

    # Find the critical potential mu_c by finding the root of the objective function.
    # The expected value based on literature is between 2 and 4.
    try:
        sol = root_scalar(objective_function, bracket=[2.0, 4.0], method='brentq')
        mu_c = sol.root
    except ValueError:
        print("Could not find a root in the given bracket. The model may behave unexpectedly or the bracket might be wrong.")
        return

    # --- Output the Final Result ---
    print("The final result is obtained by solving for 'mu' in the following differential equation:")
    # Interpolating the found value and other parameters into a string representation of the equation
    final_equation = (
        "psi''(z) + [f'(z)/f(z) - 3/z] * psi'(z) + "
        f"[{mu_c:.4f}^2*(1-z^2)^2/f(z)^2 - ({M_SQUARED})/ (z^2*f(z))] * psi(z) = 0"
    )
    print(final_equation)
    print("\nwhere f(z) is the metric function for the Einstein-Gauss-Bonnet background with:")
    print(f"Gauss-Bonnet coupling lambda_GB = {LAMBDA_GB}")
    print("\nThe numerical value found for the critical chemical potential is:")
    print(f"mu_c = {mu_c:.4f}")
    
    # Required final answer format
    print(f"\n<<<{mu_c:.4f}>>>")

# Execute the calculation
solve_critical_potential()