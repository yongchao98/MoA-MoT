import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def solve_critical_potential():
    """
    Calculates the critical chemical potential for a holographic superconductor
    in 5D Einstein-Gauss-Bonnet gravity using a shooting method.
    """
    # --- Model Parameters ---
    # The user specified the Gauss-Bonnet coupling.
    LAMBDA_GB = 0.1
    # We choose a standard value for the scalar field's mass squared,
    # corresponding to a CFT operator of dimension Delta = 3.
    M_SQ = -3.0
    DELTA_MINUS = 1.0
    DELTA_PLUS = 3.0

    # --- Background Metric Functions ---
    def f(z, lam_gb):
        """The metric function f(z) for the EGB black hole."""
        if lam_gb == 0:
            return 1.0 - z**4
        # Add a small complex part to prevent domain errors from precision issues
        sqrt_arg = 1.0 - 4.0 * lam_gb * (1.0 - z**4)
        return (1.0 - np.sqrt(sqrt_arg + 0j).real) / (2.0 * lam_gb)

    def f_prime(z, lam_gb):
        """The derivative of the metric function, f'(z)."""
        if lam_gb == 0:
            return -4.0 * z**3
        # Add a small complex part to prevent domain errors from precision issues
        sqrt_arg = 1.0 - 4.0 * lam_gb * (1.0 - z**4)
        return (4.0 * z**3) / np.sqrt(sqrt_arg + 0j).real

    # --- ODE System for the Scalar Field ---
    def odesystem(z, y, mu, lam_gb, m_sq):
        """
        Defines the system of first-order ODEs for the scalar field psi.
        The equation is psi'' + p(z)psi' + q(z)psi = 0.
        y = [psi, dpsi_dz]
        """
        psi, dpsi_dz = y
        
        f_val = f(z, lam_gb)
        fp_val = f_prime(z, lam_gb)
        
        # The background gauge field is phi(z) = mu*(1-z^2)
        phi_val = mu * (1.0 - z**2)
        
        # Coefficients of the ODE
        p_z = fp_val / f_val - 3.0 / z
        q_z = (phi_val**2 / f_val**2) + (m_sq / (z**2 * f_val))
        
        d2psi_dz2 = -p_z * dpsi_dz - q_z * psi
        return [dpsi_dz, d2psi_dz2]

    # --- Shooting Method Function ---
    def get_source_term(mu):
        """
        For a given chemical potential mu, this function solves the ODE and
        returns the source term in the asymptotic solution near the boundary.
        The critical chemical potential is the value of mu for which this source is zero.
        """
        # Integrate from near the horizon (z=1) to near the boundary (z=0)
        z_start = 1.0 - 1e-5
        z_end = 1e-6
        
        # Initial conditions at the horizon are fixed by the regularity requirement
        psi_start = 1.0  # Arbitrary normalization
        dpsi_start = -(M_SQ / 4.0) * psi_start
        
        # Solve the Initial Value Problem
        sol = solve_ivp(
            fun=odesystem,
            t_span=[z_start, z_end],
            y0=[psi_start, dpsi_start],
            args=(mu, LAMBDA_GB, M_SQ),
            method='RK45'
        )
        
        # Extract the solution at the end of the integration range (near z=0)
        psi_end = sol.y[0, -1]
        dpsi_end = sol.y[1, -1]
        
        # The asymptotic behavior near the boundary z=0 is psi = c1*z^1 + c2*z^3.
        # We need to find the value of the source term, c1.
        # From the numerical solution at z_end, we can compute c1:
        #   psi_end = c1*z_end + c2*z_end^3
        #   dpsi_end = c1 + 3*c2*z_end^2
        # Solving for c1 gives:
        c1 = (DELTA_PLUS * psi_end - z_end * dpsi_end) / ((DELTA_PLUS - DELTA_MINUS) * z_end**DELTA_MINUS)
        
        return c1

    # --- Find the Root ---
    # We search for mu_c in a physically reasonable bracket [5, 15].
    # The root of get_source_term(mu) is the critical chemical potential.
    mu_c = brentq(get_source_term, a=5.0, b=15.0)

    # --- Print the Result ---
    gauss_bonnet_coupling = 0.1
    print(f"The calculation is for a holographic model with a Gauss-Bonnet coupling of {gauss_bonnet_coupling}.")
    print("The critical chemical potential (mu_c) is found by solving the scalar field's equation of motion.")
    print("The final value is the one for which the source of the dual operator vanishes.")
    print(f"Resulting equation: mu_c = {mu_c:.4f}")

if __name__ == '__main__':
    solve_critical_potential()