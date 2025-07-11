import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def solve_critical_potential():
    """
    Calculates the critical chemical potential for a holographic superconductor
    in 5D Einstein-Gauss-Bonnet gravity.
    """
    lmbda_gb = 0.1  # Gauss-Bonnet coupling
    m2 = -3.0       # Squared mass of the scalar field (for dual operator of dimension 3)
    d = 4           # Number of boundary dimensions

    # The equation for the scalar field psi depends on the background metric.
    # The metric function f(z) for an EGB black hole (z_h=1):
    # f(z) = (1 / (2 * lmbda_gb)) * (1 - sqrt(1 - 4 * lmbda_gb * (1 - z^4)))
    def f(z, lmbda):
        # Add a small epsilon to z to avoid sqrt of negative for z > 1
        # due to floating point inaccuracies
        z = np.minimum(z, 1.0)
        return (1.0 / (2.0 * lmbda)) * (1.0 - np.sqrt(1.0 - 4.0 * lmbda * (1.0 - z**4)))

    # The derivative of the metric function, f'(z)
    def f_prime(z, lmbda):
        # Handle the limit z -> 1 where the denominator goes to zero
        if np.isclose(z, 1.0):
            # Analytically calculated limit f'(z=1) = -(d-2) = -2 in 4D if there is no r^2 infront of f(r).
            # For this form of f(z), f'(z_h) = -d/z_h = -4
            return -float(d)
        
        # Add epsilon to avoid division by zero away from the horizon
        z = np.maximum(z, 1e-9) 
        # Add small epsilon to argument of sqrt to avoid it becoming negative
        sqrt_arg = 1.0 - 4.0 * lmbda * (1.0 - z**4)
        sqrt_arg = np.maximum(sqrt_arg, 0)
        
        return (2.0 * z**3) / np.sqrt(sqrt_arg)
        
    print("In the context of a holographic superconductor in a 5D Einstein-Gauss-Bonnet background:")
    print(f"The Gauss-Bonnet coupling is set to lambda_GB = {lmbda_gb}.")
    print("The scalar field is dual to an operator of dimension Delta = 3, which corresponds to a mass-squared m^2 = -3.")

    # We need to solve the linearized equation for the scalar field psi(z):
    # psi'' + (f'/f + (d-3)/z)psi' + (mu^2*(1-z^2)^2/f^2 - m^2/(z^2*f))psi = 0
    # where d=4, f is the EGB metric function, and z_h=1.
    def ode_system(z, y, mu2, lmbda, m_squared):
        psi, psi_prime = y
        
        # Avoid division by zero at z=0
        if np.isclose(z, 0):
            return [psi_prime, 0]

        fz = f(z, lmbda)
        fpz = f_prime(z, lmbda)
        
        # Coefficients of the ODE
        coeff_psi_prime = fpz / fz + (d-3) / z
        
        # Handle the phi^2/f^2 term carefully. phi=mu(1-z^2).
        phi_term_over_f2 = mu2 * ((1-z**2)**2) / (fz**2)
        mass_term_over_f = m_squared / (z**2 * fz)
        
        coeff_psi = phi_term_over_f2 - mass_term_over_f
        
        psi_double_prime = -coeff_psi_prime * psi_prime - coeff_psi * psi
        return [psi_prime, psi_double_prime]

    # The goal is to find mu that allows a solution with the source term set to zero.
    # This is an eigenvalue problem for mu^2.
    def objective_function(mu2):
        # We shoot from near the horizon to near the boundary.
        z_start = 1.0 - 1e-5
        z_end = 1e-5

        # Initial conditions at the horizon (z=1) are fixed by regularity.
        # psi'(1) = (m^2 / f'(1)) * psi(1)
        # We can set psi(1) = 1 without loss of generality.
        psi_initial = 1.0
        psi_prime_initial = m2 / f_prime(z_start, lmbda_gb) * psi_initial
        y_initial = [psi_initial, psi_prime_initial]

        # Solve the ODE
        sol = solve_ivp(
            lambda z, y: ode_system(z, y, mu2, lmbda_gb, m2),
            [z_start, z_end],
            y_initial,
            dense_output=True,
            method='RK45'
        )
        
        # Extract the solution at the end point
        psi_final, psi_prime_final = sol.sol(z_end)
        
        # At the boundary z->0, psi behaves as psi(z) = psi_1*z + psi_3*z^3.
        # psi_1 is the source, psi_3 is the condensate VEV.
        # We want a solution with zero source (psi_1 = 0).
        # We can extract psi_1 from the numerical solution at z_end.
        # The relation is: psi_1 = (3*psi(z_end)/z_end - psi_prime(z_end))/2
        psi_1 = (3.0 * psi_final / z_end - psi_prime_final) / 2.0
        
        return psi_1

    # Find the root of the objective function to get the eigenvalue mu^2.
    # We need to provide a bracket [a,b] where the function has different signs.
    # Based on literature, we expect mu to be > 10, so mu^2 > 100.
    # The expected value for mu/T_c is ~14.7, which gives mu ~ 14.7/pi ~ 4.7
    # mu^2 is roughly 22. So let's bracket [10, 100]. Let's try to increase it based on older estimate
    # My first estimate of mu ~ 14.8 was actually mu/T * 1/pi. The correct mu = (mu/T_c) / pi. Let's adjust bracket
    mu2_lower, mu2_upper = 15, 30
    try:
        critical_mu2 = brentq(objective_function, mu2_lower, mu2_upper, xtol=1e-6, rtol=1e-6)
    except ValueError:
        print("Error: The objective function does not have different signs at the endpoints of the bracket.")
        print(f"Value at {mu2_lower}: {objective_function(mu2_lower)}")
        print(f"Value at {mu2_upper}: {objective_function(mu2_upper)}")
        print("Please adjust the bracket [mu2_lower, mu2_upper].")
        return

    critical_mu = np.sqrt(critical_mu2)
    
    # This value of mu is in units where the horizon radius z_h=1.
    # The temperature is T = 1/(pi*z_h) = 1/pi.
    # The calculated critical chemical potential mu_c corresponds to this temperature.
    # The dimensionless ratio mu_c / T_c is often the desired physical quantity.
    mu_over_T = critical_mu * np.pi
    
    print("\nTo find the critical chemical potential, we solve the following equation for the scalar field psi(z):")
    print("psi''(z) + (f'(z)/f(z) + 1/z) * psi'(z) + (mu^2*(1-z^2)^2/f(z)^2 - (-3)/(z^2*f(z))) * psi(z) = 0")
    print("where f(z) is the Einstein-Gauss-Bonnet metric function.")
    print("We solve this as an eigenvalue problem for mu^2, demanding a solution that is regular at the horizon and has zero source at the boundary.")

    print(f"\nThe numerical solver finds the lowest eigenvalue for mu^2 to be: {critical_mu2:.4f}")
    print(f"This corresponds to a critical chemical potential mu_c (in units of horizon radius z_h=1) of: sqrt({critical_mu2:.4f}) = {critical_mu:.4f}")
    
    print("\nThis value corresponds to a specific temperature T=1/pi.")
    print(f"The physically relevant dimensionless ratio is mu_c / T_c = {critical_mu:.4f} * pi = {mu_over_T:.4f}")
    print("\nHowever, the standard definition of the critical chemical potential in this context often refers to the value of this ratio.")
    print("\nThe critical chemical potential is therefore:")
    print(mu_over_T)


solve_critical_potential()