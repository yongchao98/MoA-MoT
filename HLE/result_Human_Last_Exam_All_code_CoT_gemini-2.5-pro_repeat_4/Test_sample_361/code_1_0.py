import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def solve_critical_potential():
    """
    Calculates the critical chemical potential for scalar condensation in a D3/D7 EGB model.
    """
    # Gauss-Bonnet coupling
    LGB = 0.1

    # We work in units where the horizon radius r_h = 1.
    r_h = 1.0

    # Define the metric function f(r) for the EGB black hole
    def f(r, lgb):
        # Ensure the term inside sqrt is non-negative
        sqrt_term = 1 - 4 * lgb * (1 - (r_h**4 / r**4))
        # Handle potential floating point inaccuracies at r=r_h
        sqrt_term = np.maximum(sqrt_term, 0)
        return (1 / (2 * lgb)) * (1 - np.sqrt(sqrt_term))

    # The equation of motion for the scalar fluctuation psi(r) is:
    # psi'' + P(r) * psi' + Q(r, mu) * psi = 0
    # where P and Q are coefficients derived from the model.
    def p_coeff(r, lgb):
        # P(r) = f'(r)/f(r) + 3/r
        # Derivative of f(r)
        sqrt_term = 1 - 4 * lgb * (1 - (r_h**4 / r**4))
        sqrt_term = np.maximum(sqrt_term, 0)
        # Add a small epsilon to avoid division by zero in sqrt_term if r is huge
        df_dr = (4 * lgb * r_h**4) / (r**5 * np.sqrt(sqrt_term + 1e-20))
        # Avoid division by zero at the horizon r=r_h
        f_val = f(r, lgb)
        if np.abs(f_val) < 1e-9:
            # As r -> r_h, f(r) ~ f'(r_h)*(r-r_h), so f'/f ~ 1/(r-r_h)
            f_prime_over_f = 1.0 / (r - r_h) if r != r_h else 0
        else:
            f_prime_over_f = df_dr / f_val
        return f_prime_over_f + 3.0 / r

    def q_coeff(r, mu, lgb):
        # Q(r) = (3 / (r^2*f)) + (mu^2 * (1-r_h^2/r^2)^2) / (r^4 * f^2)
        f_val = f(r, lgb)
        f_val_sq = f_val**2
        # Avoid division by zero
        if np.abs(f_val_sq) < 1e-12:
            return 0 # Value not needed right at the horizon singularity
        term1 = 3.0 / (r**2 * f_val)
        term2 = (mu**2 * (1 - r_h**2 / r**2)**2) / (r**4 * f_val_sq)
        return term1 + term2

    # Define the ODE system for solve_ivp
    # y[0] = psi, y[1] = psi'
    def ode_system(r, y, mu, lgb):
        psi, dpsi = y
        d2psi = -p_coeff(r, lgb) * dpsi - q_coeff(r, mu, lgb) * psi
        return [dpsi, d2psi]

    # Asymptotic analysis at r -> infinity
    # The metric becomes AdS with a modified radius. f(r->inf) -> c
    c = (1 / (2 * LGB)) * (1 - np.sqrt(1 - 4 * LGB))
    # The asymptotic behavior of psi is psi ~ r^(-1 +/- i*beta)
    # We need to pick one solution corresponding to the VEV (no source)
    beta = np.sqrt(3 / c - 1)
    
    # Boundary for integration
    r_max = 200.0
    r_span = [r_max, r_h]

    # Objective function for the root finder.
    # We shoot from r_max to r_h and check the horizon regularity condition.
    # The regularity condition is psi'(r_h) - (3/4)*psi(r_h) = 0.
    def objective(mu):
        # Set initial conditions at r_max for the "no source" solution.
        # This is one of the two asymptotic solutions, e.g., the sin part.
        psi_max = r_max**(-1) * np.sin(beta * np.log(r_max))
        dpsi_max = r_max**(-2) * (-np.sin(beta * np.log(r_max)) + beta * np.cos(beta * np.log(r_max)))
        y_init = [psi_max, dpsi_max]

        # Solve the ODE
        sol = solve_ivp(ode_system, r_span, y_init, args=(mu, LGB), dense_output=True, method='RK45', atol=1e-8, rtol=1e-8)

        # Get the solution at the horizon
        psi_h, dpsi_h = sol.sol(r_h)

        # Return the value of the regularity condition
        # The root of this function is the critical chemical potential
        return dpsi_h - 0.75 * psi_h

    # Find the critical chemical potential mu_c (in units of r_h)
    # Based on literature, the value for mu/r_h should be in the range [1, 5]
    try:
        mu_crit_rh_units = brentq(objective, 1.0, 5.0, xtol=1e-6)
    except ValueError:
        print("Root not found in the initial interval. The model might not have a solution or the interval is wrong.")
        return

    # The physical result is the dimensionless ratio mu_c / T
    # Temperature T is related to the horizon radius r_h
    # The factor r_h / T is pi * sqrt(1 - 4*lambda)
    rh_over_T = np.pi * np.sqrt(1 - 4 * LGB)

    # Calculate mu_c / T
    mu_c_over_T = mu_crit_rh_units * rh_over_T

    print(f"For a Gauss-Bonnet coupling of lambda_GB = {LGB}:")
    print(f"The critical chemical potential in units of the horizon radius, mu_c/r_h, is found to be: {mu_crit_rh_units:.4f}")
    print(f"The temperature is related to the horizon radius by T/r_h = 1 / (pi * sqrt(1 - 4*lambda_GB)) = 1 / {rh_over_T:.4f}")
    print("\nThe final equation for the dimensionless critical chemical potential is:")
    print(f"mu_c/T = (mu_c/r_h) * (r_h/T)")
    print(f"{mu_c_over_T:.4f} = {mu_crit_rh_units:.4f} * {rh_over_T:.4f}")
    print("\nSo, the value of the critical chemical potential, expressed as the dimensionless ratio mu_c/T, is approximately:")
    print(f"{mu_c_over_T:.4f}")
    
    # Final answer in the required format
    print(f"\n<<<{mu_c_over_T:.4f}>>>")

solve_critical_potential()