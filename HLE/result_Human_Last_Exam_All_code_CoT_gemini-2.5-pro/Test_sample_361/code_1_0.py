import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def solve_critical_potential():
    """
    Calculates the critical chemical potential for a holographic superconductor
    in 5D Einstein-Gauss-Bonnet gravity.
    """
    # --- 1. Define constants and parameters ---
    LAMBDA_GB = 0.1  # Gauss-Bonnet coupling
    M2 = -3.0        # Mass squared of the scalar field (for operator dim Delta=3)
    Q_CHARGE = 1.0   # Charge of the scalar field

    # Integration parameters
    R_H = 1.0        # Horizon radius, set to 1
    R_INF = 200.0    # Numerical infinity
    EPS = 1e-6       # Small number to avoid singularities at the horizon

    # --- 2. Define metric functions based on EGB gravity ---

    # N^2 relates bulk time to boundary time. It's also the asymptotic speed of light.
    # It must be real, which constrains LAMBDA_GB <= 1/4.
    if 1 - 4 * LAMBDA_GB < 0:
        raise ValueError("LAMBDA_GB must be <= 1/4 for a real metric solution.")
    N_SQ = 0.5 * (1 + np.sqrt(1 - 4 * LAMBDA_GB))

    def f(r, lam):
        """Metric function f(r) for the EGB black hole."""
        if r < R_H:
            return 0
        # This term can become negative if lam > 1/4, handled above.
        sqrt_term = np.sqrt(1 - 4 * lam * (1 - (R_H / r)**4))
        return (r**2 / (2 * lam)) * (1 - sqrt_term)

    def f_prime(r, lam):
        """Derivative of f(r) with respect to r."""
        if r < R_H:
            return 0
        sqrt_term = np.sqrt(1 - 4 * lam * (1 - (R_H / r)**4))
        # Handle the r=R_H case by limit, which gives 4*R_H
        if np.isclose(r, R_H):
            return 4.0 * R_H
        term1 = (r / lam) * (1 - sqrt_term)
        term2 = (4 * (R_H**4 / r**3)) / sqrt_term
        return term1 + term2

    # --- 3. Define the ODE system for the scalar field psi ---
    # The ODE is psi'' + A(r)*psi' + B(r)*psi = 0
    # where A(r) = (3/r + f'/f)
    # and B(r) = (m^2/f - mu^2*(1-R_H^2/r^2)^2 / (N_sq*f)^2)
    # We solve for y = [psi, psi']
    def ode_system(r, y, mu, lam, m2, n_sq):
        psi, dpsi = y
        
        f_r = f(r, lam)
        f_prime_r = f_prime(r, lam)

        # At large r, f(r) can be very large, causing floating point issues if not handled.
        if f_r == 0 or not np.isfinite(f_r):
            return [dpsi, 0]

        # The solution for the gauge field phi at the critical point
        phi_r = mu * (1 - R_H**2 / r**2)
        
        # Coefficients of the ODE
        A_r = (3.0 / r) + (f_prime_r / f_r)
        B_r = (m2 / f_r) - (Q_CHARGE**2 * phi_r**2) / (n_sq * f_r**2)
        
        d2psi = -A_r * dpsi - B_r * psi
        return [dpsi, d2psi]

    # --- 4. Define the shooting function ---
    # This function takes a guess for mu, solves the ODE, and returns a value
    # that should be zero when mu is the correct critical potential.
    def shoot_target_function(mu):
        """
        Solves the ODE for a given mu and returns the source term J,
        which we want to drive to zero.
        """
        # Initial conditions determined by regularity at the horizon r = R_H
        # psi'(R_H) = -m^2 * psi(R_H) / f'(R_H)
        f_prime_at_horizon = 4.0 * R_H
        psi_prime_at_horizon = -M2 / f_prime_at_horizon
        
        # Start integration slightly away from the horizon
        r_start = R_H + EPS
        
        # We can set psi(r_start) = 1 due to linearity
        psi_init = 1.0
        dpsi_init = psi_prime_at_horizon * psi_init
        
        # We need to give the initial value at r_start, not R_H
        # A simple Taylor expansion: psi(R_H+eps) ~ psi(R_H) + eps*psi'(R_H)
        y0 = [psi_init + EPS * dpsi_init, dpsi_init]

        # Solve the ODE
        sol = solve_ivp(
            fun=ode_system,
            t_span=[r_start, R_INF],
            y0=y0,
            args=(mu, LAMBDA_GB, M2, N_SQ),
            dense_output=True,
            method='RK45'
        )

        # Extract solution at the boundary r = R_INF
        psi_inf, dpsi_inf = sol.y[:, -1]

        # The asymptotic solution is psi(r) ~ J/r + O/r^3.
        # We extract the source term J. We want J=0 for condensation.
        # J = (3 * r * psi(r) + r^2 * psi'(r)) / 2
        J = (3 * R_INF * psi_inf + R_INF**2 * dpsi_inf) / 2.0
        return J

    # --- 5. Find the root to get the critical potential ---
    # We need to provide a bracket [mu_low, mu_high] where the target function changes sign.
    # Based on literature, a bracket of [1.0, 4.0] should be safe for mu_bulk.
    try:
        # We are solving for mu_bulk here.
        mu_bulk_c = brentq(shoot_target_function, a=1.0, b=4.0, xtol=1e-8)
        
        # The physical chemical potential at the boundary is mu_c = N * mu_bulk_c
        mu_c = np.sqrt(N_SQ) * mu_bulk_c
        
        # Print the final result in a descriptive sentence
        print("In a 5D Einstein-Gauss-Bonnet holographic model with:")
        print(f"  Gauss-Bonnet coupling (lambda_GB) = {LAMBDA_GB}")
        print(f"  Scalar operator dimension (Delta) = 3 (m^2 = {M2})")
        print("\nThe calculated critical chemical potential (mu_c) for condensation is:")
        print(f"mu_c = {mu_c:.6f}")
        
        return mu_c

    except ValueError:
        print("Failed to find a root. The function may not cross zero in the given interval.")
        print("Please try adjusting the bracket [a, b] for the root finder.")
        return None

# Execute the calculation and store the final answer
critical_potential = solve_critical_potential()
if critical_potential is not None:
    print(f"\n<<< {critical_potential:.6f} >>>")