import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar

def solve_critical_potential():
    """
    Calculates the critical chemical potential for scalar condensation in a
    holographic D3/D7 model with an Einstein-Gauss-Bonnet background.
    """
    # Gauss-Bonnet coupling constant
    LAM_GB = 0.1

    # Define the metric function f(z) and its derivative f'(z)
    # A numerically stable form of f(z) is used.
    def f(z, lam_gb):
        if z == 0:
            return np.inf
        if z == 1:
            return 0
        sqrt_term = np.sqrt(1 - 4 * lam_gb * (1 - z**4))
        # This form is stable for small lam_gb
        return 2 * (1 - z**4) / (z**2 * (1 + sqrt_term))

    def f_prime(z, lam_gb):
        if z == 0:
            return np.inf
        # The derivative can be calculated analytically.
        # f'(z) = d/dz [ (1 - sqrt(1 - 4*lam_gb*(1-z^4))) / (2*lam_gb*z^2) ]
        sqrt_term = np.sqrt(1 - 4 * lam_gb * (1 - z**4))
        numerator = -16 * lam_gb**2 * z**5 / sqrt_term - (1 - sqrt_term) * 4 * lam_gb * z
        denominator = (2 * lam_gb * z**2)**2
        return numerator / denominator

    # Define the system of first-order ODEs for the solver
    # y[0] = Psi, y[1] = Psi'
    def ode_system(z, y, mu, lam_gb):
        Psi, dPsi_dz = y[0], y[1]
        
        fz = f(z, lam_gb)
        fpz = f_prime(z, lam_gb)

        # Coefficients of the ODE: Psi'' + P(z)Psi' + Q(z)Psi = 0
        P_z = fpz / fz - 3 / z
        
        # Mass of the scalar field is m^2 = -3
        m_squared = -3
        Q_z = (m_squared / (z**2 * fz)) + (mu**2 * (1 - z**2)**2) / (z**2 * fz**2)
        
        d2Psi_dz2 = -P_z * dPsi_dz - Q_z * Psi
        
        return [dPsi_dz, d2Psi_dz2]

    # The "shooting" function. It integrates the ODE for a given mu
    # and returns the value of Psi at the boundary.
    # We want to find mu such that this function is zero.
    def shoot_objective(mu):
        # Start integration slightly away from the horizon singularity at z=1
        z_start = 1 - 1e-5
        # End integration close to the boundary at z=0
        z_end = 1e-5

        # Set initial conditions using the horizon regularity condition
        # Psi'(1) = (m_squared / f'(1)) * Psi(1). Here m_squared = -3.
        # We can set Psi(1) = 1 due to linearity.
        f_prime_at_1 = 4 / np.sqrt(1 - 4 * LAM_GB)
        
        psi_0 = 1.0
        dpsi_dz_0 = (m_squared / f_prime_at_1) * psi_0

        # Solve the initial value problem
        sol = solve_ivp(
            fun=ode_system,
            t_span=[z_start, z_end],
            y0=[psi_0, dpsi_dz_0],
            args=(mu, LAM_GB),
            dense_output=True,
            method='RK45'
        )
        
        # Return the value of Psi at the endpoint.
        # The root finder will find mu that makes this zero.
        return sol.y[0, -1]

    # Use a root-finding algorithm to find the critical chemical potential mu_c
    # We need to provide a bracket [a, b] where shoot_objective(a) and shoot_objective(b) have opposite signs.
    # Based on known physics, the value should be around 1.
    try:
        result = root_scalar(shoot_objective, bracket=[0.5, 1.5], method='brentq')
        mu_c = result.root
        
        print(f"Calculation for Gauss-Bonnet coupling lambda_GB = {LAM_GB}")
        print(f"The scalar operator has dimension Delta=3, corresponding to mass^2 = -3.")
        print("\nThe final equation being solved numerically is of the form:")
        print("Psi''(z) + P(z) * Psi'(z) + Q(z, mu) * Psi(z) = 0")
        print(f"The critical chemical potential is found to be:")
        print(f"mu_c = {mu_c:.8f}")
        
        # Final answer in the required format
        print(f"\n<<<{mu_c:.8f}>>>")

    except ValueError as e:
        print(f"Error during root finding: {e}")
        print("The initial bracket might not contain a root. Please try adjusting the bracket [0.5, 1.5].")

if __name__ == '__main__':
    solve_critical_potential()