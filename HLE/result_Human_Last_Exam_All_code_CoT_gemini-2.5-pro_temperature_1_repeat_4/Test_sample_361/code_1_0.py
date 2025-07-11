import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def solve_critical_potential():
    """
    Calculates the critical chemical potential for scalar condensation in a 5D
    Einstein-Gauss-Bonnet holographic model.
    """

    # --- 1. Model Parameters ---
    # The user specifies the Gauss-Bonnet coupling.
    lambda_gb = 0.1
    # In the D3/D7 model, the operator dual to the condensing scalar is the quark
    # bilinear <bar{q}q>, which has conformal dimension Delta = 3.
    # The mass-squared of the dual scalar field in the bulk is related to Delta by
    # m^2 = Delta * (Delta - d), where d=4 is the boundary spacetime dimension.
    m2 = 3 * (3 - 4)

    # --- 2. Define the metric function and its derivative ---
    # The metric is ds^2 = 1/z^2 * (-f(z)dt^2 + dx_i^2 + 1/f(z)dz^2)
    # with z_h=1 (horizon radius). z=0 is the boundary.
    # We add a small epsilon to the term under the square root to prevent
    # numerical errors right at z=1, although the integration will start just shy of 1.
    epsilon = 1e-12
    def f(z, lgb):
        """Metric function f(z) for EGB black hole."""
        sqrt_arg = 1 - 4 * lgb * (1 - z**4)
        return (1 / (2 * lgb)) * (1 - np.sqrt(sqrt_arg + epsilon))

    def f_prime(z, lgb):
        """Derivative of f(z) with respect to z."""
        sqrt_arg = 1 - 4 * lgb * (1 - z**4)
        return -(4 * z**3) / np.sqrt(sqrt_arg + epsilon)

    # --- 3. Define the system of ODEs ---
    # We are solving the equation for the scalar field psi at the critical point.
    # The equation is: psi'' + (f'/f - 3/z)psi' + (mu^2(1-z^2)^2/f^2 - m^2/(z^2*f))psi = 0
    # We convert this to a system of two first-order ODEs:
    # y[0] = psi, y[1] = chi = psi'
    def odes(z, y, mu, lgb, m2_val):
        """The system of first-order ODEs for psi and psi'."""
        psi, chi = y
        
        fz = f(z, lgb)
        fpz = f_prime(z, lgb)

        coeff_chi = fpz / fz - 3.0 / z
        coeff_psi = (mu**2 * (1 - z**2)**2) / fz**2 - m2_val / (z**2 * fz)

        d_psi_dz = chi
        d_chi_dz = -coeff_chi * chi - coeff_psi * psi

        return [d_psi_dz, d_chi_dz]

    # --- 4. Define the objective function for the shooting method ---
    def objective_function(mu, lgb, m2_val):
        """
        Solves the ODE for a given mu and returns the value of psi at the boundary.
        The goal is to find mu such that this function returns 0.
        """
        z_start = 1.0 - 1e-6
        z_end = 1e-6

        # Initial conditions from regularity at the horizon z=1.
        psi0 = 1.0
        chi0 = -m2_val / 4.0
        y0 = [psi0, chi0]

        # Solve the ODE system
        sol = solve_ivp(
            fun=lambda z, y: odes(z, y, mu, lgb, m2_val),
            t_span=[z_start, z_end],
            y0=y0,
            method='Radau',  # A good method for stiff ODEs
            dense_output=True
        )

        # Value of psi at the end of the integration (near z=0).
        # Normalizability requires this to be zero.
        psi_at_boundary = sol.sol(z_end)[0]
        return psi_at_boundary

    # --- 5. Find the critical chemical potential using a root finder ---
    # We use brentq, which requires a bracket [a, b] where f(a) and f(b) have opposite signs.
    # From literature, for lambda_gb=0.1, mu_c should be around 3.65.
    lower_bound = 3.0
    upper_bound = 4.0

    try:
        critical_mu = brentq(
            f=objective_function,
            a=lower_bound,
            b=upper_bound,
            args=(lambda_gb, m2)
        )
        # The prompt asks to output each number in the final equation.
        # The final "equation" can be stated as:
        # mu_c(lambda_gb, m^2) = result
        print(f"For a Gauss-Bonnet coupling of {lambda_gb} and a scalar mass-squared of {m2},")
        print(f"the critical chemical potential is {critical_mu:.4f}")

    except ValueError:
        print("The root-finding bracket is not valid.")
        print("The objective function may not have opposite signs at the bounds.")
        print(f"Value at lower bound ({lower_bound}): {objective_function(lower_bound, lambda_gb, m2)}")
        print(f"Value at upper bound ({upper_bound}): {objective_function(upper_bound, lambda_gb, m2)}")
        
solve_critical_potential()