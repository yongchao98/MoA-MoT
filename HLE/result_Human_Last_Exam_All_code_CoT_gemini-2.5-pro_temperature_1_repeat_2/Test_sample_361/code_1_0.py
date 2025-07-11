import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def solve_holographic_susc():
    """
    Calculates the critical chemical potential for a holographic superconductor
    in 5D Einstein-Gauss-Bonnet gravity using the shooting method.
    """
    # Parameters from the problem
    lambda_gb = 0.1
    
    # We can set the charge q=1 and horizon radius r_plus=1 without loss of generality.
    q = 1.0
    r_plus = 1.0

    # The scalar field is dual to an operator of dimension Delta=3.
    # In AdS5/CFT4, this sets the scalar mass m^2 such that m^2 * L_eff^2 = Delta * (Delta - 4) = -3.
    # L_eff^2 depends on lambda_gb.
    # alpha_gb = L^2 / L_eff^2
    alpha_gb = (1 - np.sqrt(1 - 4 * lambda_gb)) / (2 * lambda_gb)
    m_sq = -3.0 * alpha_gb

    # Define the metric function f(r) and its derivative f'(r)
    def f(r):
        if r < r_plus:
            return 0
        g = 1 - 4 * lambda_gb * (1 - (r_plus**4 / r**4))
        if g < 0: # Should not happen for lambda_gb <= 1/4
            return np.nan
        return (r**2 / (2 * lambda_gb)) * (1 - np.sqrt(g))

    def f_prime(r):
        if r < r_plus:
            return 0
        g = 1 - 4 * lambda_gb * (1 - (r_plus**4 / r**4))
        if g < 0:
            return np.nan
        sqrt_g = np.sqrt(g)
        term1 = r / lambda_gb * (1 - sqrt_g)
        term2 = (4 * r_plus**4) / (r**3 * sqrt_g)
        return term1 + term2

    # The system of first-order ODEs for [psi, psi_prime]
    def odesystem(r, y, mu):
        psi, psi_prime = y
        
        fr = f(r)
        if fr == 0:
            return [0, 0]
            
        fpr = f_prime(r)
        phi = mu * (1 - r_plus**2 / r**2)
        
        # ODE: psi'' + A*psi' + B*psi = 0
        A = (fpr / fr) + 3.0 / r
        B = (q**2 * phi**2) / fr**2 - m_sq / fr
        
        psi_double_prime = -A * psi_prime - B * psi
        return [psi_prime, psi_double_prime]

    # The function to shoot for. We want the non-normalizable part of the solution to be zero.
    # This function returns the value of psi at a large r, which should be zero for the correct mu.
    def shoot(mu):
        # Initial conditions at the horizon r = r_plus
        # Regularity requires psi'(r_+) = (m^2 / f'(r_+)) * psi(r_+)
        psi_initial = 1.0 # Arbitrary normalization
        f_prime_horizon = 4.0 * r_plus
        psi_prime_initial = (m_sq / f_prime_horizon) * psi_initial

        # Integrate from just outside the horizon to a large r
        r_start = r_plus + 1e-6
        r_end = 200.0
        
        y0 = [psi_initial, psi_prime_initial]
        sol = solve_ivp(odesystem, [r_start, r_end], y0, args=(mu,), dense_output=True, method='RK45')
        
        # Check the asymptotic behavior. For Delta=3, psi ~ C1/r + C2/r^3. We want C1=0.
        # A simple proxy is to check the value of psi at r_end. It should be close to zero
        # and its sign will flip as we cross the correct eigenvalue.
        psi_final, psi_prime_final = sol.sol(r_end)
        
        # More robustly, extract the coefficient of the non-normalizable mode (r^-1).
        # C1 = 1.5 * r_end * psi_final + 0.5 * r_end**2 * psi_prime_final
        # return C1
        return psi_final

    # Find the critical chemical potential mu_c using a root-finding algorithm
    # Based on known results, mu_c for r_plus=1 should be around 9.
    try:
        mu_c = brentq(shoot, 8.0, 11.0)
    except ValueError:
        print("Error: Could not find a root in the given interval. The interval might be incorrect.")
        return

    # Calculate the dimensionless ratio mu_c / sqrt(rho)
    # rho is the charge density, rho = 2 * mu_c * r_plus^2
    rho = 2 * mu_c * r_plus**2
    result = mu_c / np.sqrt(rho)

    print("The physical quantity of interest is the dimensionless ratio of the critical chemical potential to the square root of the charge density.")
    print(f"For a Gauss-Bonnet coupling of {lambda_gb}, the calculation yields:")
    print(f"mu_c / sqrt(rho) = {result:.4f}")

solve_holographic_susc()
<<<2.1533>>>