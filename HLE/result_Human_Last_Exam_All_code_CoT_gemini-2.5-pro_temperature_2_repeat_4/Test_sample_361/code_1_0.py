import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import warnings

# Suppress warnings that can occur during numerical integration for cleaner output
warnings.filterwarnings("ignore", category=np.ComplexWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

def calculate_critical_potential():
    """
    This function calculates the critical chemical potential for scalar condensation
    in a 5D Einstein-Gauss-Bonnet holographic model using a numerical shooting method.
    """

    # --- Problem Parameters ---
    # Gauss-Bonnet coupling
    lambda_gb = 0.1
    # Mass squared of the scalar field. Corresponds to a Delta=3 operator (e.g., <q-bar q>)
    # via the holographic relation m^2 = Delta(Delta-d) for d=4.
    m_sq = -3.0
    # We work in units where the black hole horizon is at zh = 1.0, which sets the energy scale.
    zh = 1.0

    # --- Define the metric function f(z) and its derivative f'(z) ---
    def f(z, lgb):
        # Handle the lgb -> 0 limit, which is the AdS-Schwarzschild case
        if lgb == 0:
            return 1 - z**4
        # For small z, use Taylor expansion f(z) ~ 1 - z^4 to avoid numerical errors
        if z < 1e-9:
            return 1.0 - z**4
        # General case for f(z) in EGB gravity
        term_in_sqrt = 1 - 4 * lgb * (1 - z**4)
        if term_in_sqrt < 0:
            return np.nan
        return (1 - np.sqrt(term_in_sqrt)) / (2 * lgb)

    def f_prime(z, lgb):
        # Handle the lgb -> 0 limit
        if lgb == 0:
            return -4 * z**3
        # Small z behavior is proportional to z^3
        if z < 1e-9:
            return -4 * z**3 * (1 - 2*lgb) # First order correction
        # General case for f'(z)
        term_in_sqrt = 1 - 4 * lgb * (1 - z**4)
        if term_in_sqrt <= 0: # The derivative is undefined if sqrt is zero
            return np.nan
        return (4 * z**3) / np.sqrt(term_in_sqrt)

    # --- Define the ODE system for the scalar field psi(z) ---
    # The system is y' = [psi', psi''], with y[0] = psi, y[1] = psi'
    def ode_system(z, y, mu_val, lgb, m_sq_val):
        psi, psi_prime = y

        if z < 1e-9: # Start integration just after z=0
            return [y[1], 0.0]

        fz = f(z, lgb)
        fz_p = f_prime(z, lgb)

        if np.isnan(fz) or np.isnan(fz_p) or fz == 0:
            return [np.nan, np.nan] # Signal solver to stop if out of domain

        # ODE coefficients P(z) and Q(z) from the main equation
        P_z = fz_p / fz - 3 / z
        # Note: (1-z/zh)^2 becomes (1-z)^2 since we set zh=1
        Q_z = (mu_val**2 * (1 - z)**2) / (fz**2) - m_sq_val / (z**2 * fz)

        d_psi_prime_dz = -P_z * psi_prime - Q_z * psi
        return [psi_prime, d_psi_prime_dz]

    # --- Function for the shooting method ---
    # This function returns a "residue" which is zero for the correct mu_c.
    def shoot(mu_val, lgb=lambda_gb, m_sq_val=m_sq, end_point=zh):
        # Start integration infinitesimally close to the boundary z=0
        epsilon = 1e-7
        # Conformal dimension of the dual operator, Delta+. For m^2=-3, Delta+ = 3.
        Delta_plus = 2 + np.sqrt(4 + m_sq_val)

        # Set initial conditions based on the asymptotic behavior psi ~ z^Delta+
        y0 = [epsilon**Delta_plus, Delta_plus * epsilon**(Delta_plus - 1)]

        # Numerically integrate the ODE
        sol = solve_ivp(
            lambda z, y: ode_system(z, y, mu_val, lgb, m_sq_val),
            [epsilon, end_point], y0, dense_output=True,
            method='RK45', rtol=1e-8, atol=1e-10
        )

        if sol.status != 0:
            return np.nan # Integration failed

        # Evaluate the solution near the horizon to avoid the z=zh singularity
        z_eval = end_point * (1 - 1e-6)
        psi_at_eval, psi_prime_at_eval = sol.sol(z_eval)

        # Regularity condition at the horizon: psi'(zh) - C * psi(zh) = 0
        f_prime_at_zh = -4.0 / (end_point * np.sqrt(1 - 4 * lgb))
        C = m_sq_val / f_prime_at_zh
        residue = psi_prime_at_eval - C * psi_at_eval
        return residue

    # --- Find mu_c using a root-finding algorithm ---
    # We search for mu_c in a plausible range. Based on literature, for lgb=0.1, it's around 4.5.
    mu_min, mu_max = 4.0, 5.0
    try:
        # brentq is an efficient and robust root-finding algorithm
        mu_c = brentq(shoot, mu_min, mu_max, xtol=1e-7, rtol=1e-7)
    except (ValueError, RuntimeError):
        print("Root could not be found in the bracket [4.0, 5.0]. The solution may lie outside this range.")
        return None
    
    # --- Print Results ---
    print("--- Holographic Condensation Calculation ---")
    print("This script solves for the critical chemical potential (mu_c) by finding the eigenvalue of the scalar field equation of motion.\n")
    print("The solved equation is: d^2(psi)/dz^2 + P(z) * d(psi)/dz + Q(z, mu_c) * psi = 0")
    print("\nWith functions defined as:")
    print(f"  P(z) = f'(z)/f(z) - 3/z")
    print(f"  Q(z, mu_c) = (mu_c^2 * (1-z)^2)/f(z)^2 - m^2/(z^2 * f(z))")
    print(f"  f(z) = (1 - sqrt(1 - 4*lambda_GB*(1 - z^4))) / (2*lambda_GB)")
    print("\nNumerical values used in the equation:")
    print(f"  Gauss-Bonnet coupling lambda_GB = {lambda_gb}")
    print(f"  Scalar mass squared m^2 = {m_sq}")
    print(f"  Horizon radius zh = {zh} (This sets the overall energy scale)")

    print("\n--- Result ---")
    print("The value of the critical chemical potential for condensation is:")
    # We explicitly print the numbers used in the equation one last time for clarity.
    print(f"mu_c (for lambda_GB={lambda_gb}, m^2={m_sq}) = {mu_c:.6f}")
    
    return mu_c

if __name__ == '__main__':
    final_answer = calculate_critical_potential()
    if final_answer:
        # This part formats the final numeric answer for the platform
        print(f"\n<<< {final_answer:.6f} >>>")
