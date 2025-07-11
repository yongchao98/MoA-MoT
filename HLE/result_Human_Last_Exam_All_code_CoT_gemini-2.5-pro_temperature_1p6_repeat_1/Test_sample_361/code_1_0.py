import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def solve_critical_potential():
    """
    Calculates the critical chemical potential for a holographic superconductor
    in 5D Einstein-Gauss-Bonnet gravity using a shooting method.
    """
    # 1. Define physical and numerical parameters
    L_GB = 0.1  # Gauss-Bonnet coupling
    M_SQ = -3.0 # Scalar field mass squared, for a Delta=3 operator
    Q_CHARGE = 1.0    # Scalar field charge

    R_START = 1.0 + 1e-7 # Start integration just outside the horizon
    R_END = 2000.0       # A large radius representing infinity

    # 2. Define metric function f(r) and its derivative f'(r) for the EGB background
    def f(r, l_gb):
        if r <= 1.0:
            return 0.0
        # Protect against float precision errors making the argument negative
        term_in_sqrt = max(0, 1.0 - 4.0 * l_gb * (1.0 - 1.0 / r**4))
        return (r**2 / (2.0 * l_gb)) * (1.0 - np.sqrt(term_in_sqrt))

    def f_prime(r, l_gb):
        if r == 1.0:
            # The value is known from a series expansion around r=1
            return 4.0
        term_in_sqrt = max(0, 1.0 - 4.0 * l_gb * (1.0 - 1.0 / r**4))
        # This case is to prevent division by zero, should not be hit in normal integration
        if term_in_sqrt == 0:
            return np.inf
        sqrt_term = np.sqrt(term_in_sqrt)
        term1 = r / l_gb * (1.0 - sqrt_term)
        term2 = (r**2 / (2.0 * l_gb)) * (8.0 * l_gb / (r**5 * sqrt_term))
        return term1 + term2

    # 3. Define the system of ODEs for psi(r)
    def ode_system(r, y, mu, l_gb, m_sq, q):
        psi, dpsi = y
        
        fr = f(r, l_gb)
        # Avoid division by zero, this check is for safety
        if fr == 0:
            return [dpsi, 0.0]
        
        f_pr = f_prime(r, l_gb)
        # The background electric potential
        phi = mu * (1.0 - 1.0 / r**2)
        
        # The linearized ODE for the scalar field psi:
        # psi'' + (f'/f + 3/r)psi' + (q^2*phi^2/f^2 - m^2/f)psi = 0
        d2psi = -(f_pr / fr + 3.0 / r) * dpsi - (q**2 * phi**2 / fr**2 - m_sq / fr) * psi
        return [dpsi, d2psi]

    # 4. Define the shooting function
    def shooting_function(mu):
        # We need to calculate the fall-off exponents at infinity to extract the source term
        alpha_sq = (1 - np.sqrt(1 - 4 * L_GB)) / (2 * L_GB)
        m_eff_sq = M_SQ / alpha_sq
        sqrt_term_p = np.sqrt(4 + m_eff_sq)
        p1 = -2 + sqrt_term_p  # Corresponds to the non-normalizable mode (source)
        p2 = -2 - sqrt_term_p  # Corresponds to the normalizable mode (VEV)

        # Initial conditions from regularity at the horizon r=1
        psi0 = 1.0
        dpsi0 = M_SQ / f_prime(1.0, L_GB)
        y0 = [psi0, dpsi0]
        
        # Numerically solve the ODE
        sol = solve_ivp(
            ode_system,
            [R_START, R_END],
            y0,
            args=(mu, L_GB, M_SQ, Q_CHARGE),
            method='LSODA',
            atol=1e-10,
            rtol=1e-10
        )
        
        # Extract the source term coefficient from the solution at R_END
        r_inf = sol.t[-1]
        psi_inf, dpsi_inf = sol.y[:, -1]
        
        # The solution at infinity is psi(r) ~ c1 * r^p1 + c2 * r^p2. We want to find mu such that c1=0.
        # We can find c1 from a linear combination of psi and psi'
        # c1 = r^(1-p1) * (r*psi' - p2*psi) / (p1 - p2)
        c1 = (r_inf**(1 - p1)) * (r_inf * dpsi_inf - p2 * psi_inf) / (p1 - p2)
        return c1

    # 5. Use a root-finder to get the critical chemical potential
    # Based on similar problems in the literature, we search in the interval [4, 7]
    try:
        mu_c = brentq(shooting_function, a=4.0, b=7.0, xtol=1e-7)
        # Final Output
        print("Final Equation Parameters:")
        print(f"The critical chemical potential (mu_c) is calculated for the following system:")
        print(f"  Gauss-Bonnet coupling (lambda_GB) = {L_GB}")
        print(f"  Scalar field mass squared (m^2) = {M_SQ}")
        print(f"Resulting critical chemical potential:")
        print(f"  mu_c = {mu_c:.6f}")
    except ValueError:
        print("Error: Could not find the critical potential. The shooting function might not cross zero in the initial bracket.")
        
solve_critical_potential()