import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def solve_critical_potential():
    """
    Calculates the critical chemical potential for a holographic superconductor
    in 5D Einstein-Gauss-Bonnet gravity using a shooting method.
    """
    # 1. Define the physical parameters of the model.
    # The user specified the Gauss-Bonnet coupling.
    LAMBDA_GB = 0.1
    # We use the standard mass for a scalar field dual to a dimension 3 operator.
    M2L2 = -3

    # Define numerical precision parameters
    Z_BOUNDARY = 1e-6  # How close we integrate towards the boundary z=0
    Z_HORIZON_OFFSET = 1e-6 # How close we start integrating from the horizon z=1
    
    # Define a search interval for the root finder. Based on known results in the literature,
    # the value should be around 8 for these parameters.
    MU_MIN = 7.0
    MU_MAX = 9.0

    # 2. Define helper functions for the EGB metric.
    def egb_f(z, lambda_gb):
        """Metric function f(z) for 5D Einstein-Gauss-Bonnet gravity."""
        # Use Taylor expansion near z=1 to avoid numerical instability
        if 1 - z**4 < 1e-12:
            return 4 * (1 - z)
        
        inner_sqrt = 1 - 4 * lambda_gb * (1 - z**4)
        return (1 - np.sqrt(inner_sqrt)) / (2 * lambda_gb)

    def egb_f_prime(z, lambda_gb):
        """Derivative of the metric function f(z)."""
        # The derivative at z=1 is exactly -4 for any valid lambda_gb.
        if 1 - z**4 < 1e-12:
            return -4.0
        
        inner_sqrt = 1 - 4 * lambda_gb * (1 - z**4)
        return 4 * z**3 / np.sqrt(inner_sqrt)

    # 3. Define the ODE system for the scalar field psi.
    def scalar_ode_system(z, y, mu, lambda_gb, m2l2):
        """
        Defines the system of first-order ODEs for the scalar field psi.
        The state vector is y = [psi, dpsi/dz].
        The function returns [dpsi/dz, d2psi/dz2].
        """
        psi, dpsi = y
        
        f_val = egb_f(z, lambda_gb)
        fp_val = egb_f_prime(z, lambda_gb)

        # This is the linearized equation of motion for the scalar field psi at the critical point.
        # It is of the form: psi'' + P(z)psi' + Q(z, mu)psi = 0
        d2psi = -(fp_val / f_val - 3 / z) * dpsi - \
                (mu**2 * (1 - z**2)**2 / f_val**2 - m2l2 / (z**2 * f_val)) * psi
                
        return [dpsi, d2psi]

    # 4. Define the shooting function.
    def find_source_term(mu):
        """
        This is the core of the shooting method. It solves the ODE for a given mu
        and returns the source term `A` from the boundary expansion psi ~ A*z + B*z^3.
        We are looking for the value of mu that makes A=0.
        """
        # The exponent for the regular solution at the horizon is nu = sqrt(-m^2 L^2 / 4).
        horizon_exponent = np.sqrt(-M2L2 / 4.0)

        # Set initial conditions near the horizon z=1, ensuring a regular solution.
        z_start = 1 - Z_HORIZON_OFFSET
        # The solution psi(z) behaves like (1-z)^nu near the horizon.
        psi_start = Z_HORIZON_OFFSET**horizon_exponent
        # Its derivative dpsi/dz behaves like -nu * (1-z)^(nu-1).
        dpsi_start = -horizon_exponent * Z_HORIZON_OFFSET**(horizon_exponent - 1)
        
        y_start = [psi_start, dpsi_start]
        
        # Solve the ODE from near the horizon to near the boundary.
        sol = solve_ivp(
            scalar_ode_system,
            [z_start, Z_BOUNDARY],
            y_start,
            args=(mu, LAMBDA_GB, M2L2),
            method='RK45',
            dense_output=True,
            atol=1e-10, rtol=1e-10
        )
        
        # Extract the solution at the final point near the boundary.
        z_end = Z_BOUNDARY
        psi_end, dpsi_end = sol.sol(z_end)
        
        # From the boundary behavior psi(z) = A*z + B*z^3, we can extract the source term A.
        source_term_A = (3 * psi_end / z_end - dpsi_end) / 2.0
        
        return source_term_A

    # 5. Find the critical potential using a root finder.
    try:
        # It is good practice to check that the function crosses zero in the interval.
        if np.sign(find_source_term(MU_MIN)) == np.sign(find_source_term(MU_MAX)):
             print("Error: Search interval is not valid. The function does not cross zero.")
             return

        # `brentq` is a robust and efficient root-finding algorithm.
        mu_c = brentq(find_source_term, MU_MIN, MU_MAX, xtol=1e-7)
        
        # 6. Output the final results.
        print("This script calculates the critical chemical potential (mu_c) for a holographic superconductor model.")
        print("The model is based on 5D Einstein-Gauss-Bonnet gravity with a scalar field dual to a dimension 3 operator.")
        
        print("\nThe parameters used in the final equation are:")
        print(f"  - Gauss-Bonnet coupling (lambda_GB): {LAMBDA_GB}")
        print(f"  - Scalar field mass squared (m^2 L^2): {M2L2}")
        
        print(f"\nThe calculated critical chemical potential is:")
        print(f"mu_c = {mu_c:.5f}")

    except Exception as e:
        print(f"An error occurred during the calculation: {e}")


if __name__ == '__main__':
    solve_critical_potential()
