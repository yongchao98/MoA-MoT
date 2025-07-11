import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar

def solve_critical_potential():
    """
    Calculates the critical chemical potential for a holographic superconductor
    in 5D Gauss-Bonnet gravity using a shooting method.
    """
    # --- Model Parameters ---
    # Gauss-Bonnet coupling
    LAM_GB = 0.1
    # Mass squared of the scalar field (for dual operator with Delta=3)
    M2 = -3.0
    # Charge of the scalar field (can be set to 1 without loss of generality)
    Q = 1.0
    # Horizon radius (sets the energy scale, r_h=1)
    RH = 1.0

    # --- EGB metric functions ---
    def g(r, lam_gb):
        """ The EGB black hole metric function f(r), denoted here as g(r). """
        # The argument of the square root must be non-negative.
        # This gives the physical constraint lambda_GB <= 1/4.
        arg_sqrt = 1 - 4 * lam_gb * (1 - (RH**4 / r**4))
        if arg_sqrt < 0:
            raise ValueError("Invalid value of lambda_GB or r, sqrt argument is negative.")
        return (r**2 / (2 * lam_gb)) * (1 - np.sqrt(arg_sqrt))

    def g_prime(r, lam_gb):
        """ The derivative of g(r) with respect to r. """
        # This form is numerically stable near the horizon.
        arg_sqrt = 1 - 4 * lam_gb * (1 - (RH**4 / r**4))
        if arg_sqrt < 0:
            raise ValueError("Invalid value of lambda_GB or r, sqrt argument is negative.")
        return (2 * g(r, lam_gb) / r) + (4 * RH**4 / (r**3 * np.sqrt(arg_sqrt)))

    # --- ODE System ---
    def ode_system(r, y, mu, lam_gb):
        """
        Defines the system of first-order ODEs for the scalar field phi.
        y = [phi, phi_prime]
        """
        phi, phi_prime = y
        
        g_r = g(r, lam_gb)
        g_prime_r = g_prime(r, lam_gb)
        
        # A_t(r) profile in the normal (non-condensed) phase
        A_t = mu * (1 - RH**2 / r**2)
        
        # The ODE is: phi'' + P(r)phi' + Q(r)phi = 0
        P_r = g_prime_r / g_r + 3.0 / r
        Q_r = (Q**2 * A_t**2) / (g_r**2) - M2 / g_r
        
        d_phi_dr = phi_prime
        d_phi_prime_dr = -P_r * phi_prime - Q_r * phi
        
        return [d_phi_dr, d_phi_prime_dr]

    # --- Shooting Method Objective Function ---
    def objective_function(mu, lam_gb, r_max, r_final):
        """
        This function represents the condition that must be met at the horizon.
        We shoot from the boundary (r_max) and check the solution near the horizon (r_final).
        The root of this function is the desired critical chemical potential mu_c.
        """
        # Boundary conditions at r_max >> RH for a condensate (phi ~ 1/r^3)
        # We set the overall normalization constant to 1.
        phi_inf = 1.0 / r_max**3
        phi_prime_inf = -3.0 / r_max**4
        y_inf = [phi_inf, phi_prime_inf]
        
        # Integrate the ODE from r_max down to r_final
        sol = solve_ivp(
            fun=lambda r, y: ode_system(r, y, mu, lam_gb),
            t_span=[r_max, r_final],
            y0=y_inf,
            dense_output=True,
            method='RK45',
            atol=1e-8,
            rtol=1e-8
        )
        
        # Extract the solution at the final point near the horizon
        phi_final, phi_prime_final = sol.sol(r_final)
        
        # Horizon regularity condition: phi'(RH) - (M2/4)*phi(RH) = 0
        # For M2 = -3, this is phi'(RH) + (3/4)*phi(RH) = 0
        error = phi_prime_final + (3.0 / 4.0) * phi_final
        
        return error

    # --- Find the Root ---
    # Integration parameters
    r_max = 1000.0  # A large value for the "boundary"
    r_final = RH + 1e-4  # A point very close to the horizon
    
    # Bracket for the root search, based on known results in the literature
    mu_bracket = [0.5, 2.0]
    
    # Use a root-finding algorithm to find mu_c
    solution = root_scalar(
        lambda mu: objective_function(mu, LAM_GB, r_max, r_final),
        bracket=mu_bracket,
        method='brentq'
    )
    mu_c = solution.root

    # --- Output the Result ---
    # The prompt asks to output each number in the "final equation".
    # The final result is a function of the input parameters: mu_c(lambda_GB, m^2).
    # We print all these numbers.
    print(f"For the holographic model defined by the parameters:")
    print(f"Gauss-Bonnet coupling, lambda_GB = {LAM_GB}")
    print(f"Scalar field mass squared, m^2 = {M2}")
    print("\nThe resulting critical chemical potential is:")
    print(f"mu_c = {mu_c:.5f}")

    return mu_c

if __name__ == '__main__':
    critical_potential = solve_critical_potential()
    # The final answer tag required by the prompt
    print(f"\n<<<{critical_potential:.5f}>>>")