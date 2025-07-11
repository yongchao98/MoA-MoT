import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar

def calculate_critical_potential():
    """
    Calculates the critical chemical potential for scalar condensation in a 
    holographic D3/D7 model with an Einstein-Gauss-Bonnet background.
    """
    
    # Key parameters defining the physics problem
    # The Gauss-Bonnet coupling as specified in the problem
    lambda_gb = 0.1
    # Effective mass squared of the scalar field in this model
    m_sq = -1.0
    # Squared charge of the scalar field
    q_sq = (0.5)**2

    # --- Define the metric function f(z) and its derivative f'(z) ---
    def f(z, lgb):
        """ Metric function f(z) for the EGB black hole. """
        # The term inside the square root must be non-negative.
        # For λ_GB=0.1 and z in [0,1], this is always satisfied.
        sqrt_term = np.sqrt(1 - 4 * lgb * (1 - z**4))
        return (1 - sqrt_term) / (2 * lgb)

    def fp(z, lgb):
        """ Derivative of the metric function, f'(z). """
        if z == 0:
            return 0.0
        sqrt_term = np.sqrt(1 - 4 * lgb * (1 - z**4))
        return (4 * z**3) / sqrt_term

    # --- Define the Ordinary Differential Equation (ODE) system ---
    def ode_system(z, y, mu, lgb, m_sq_val, q_sq_val):
        """
        Defines the system of first-order ODEs for the scalar field χ.
        y[0] = χ(z), y[1] = χ'(z).
        The function returns [χ', χ''].
        """
        chi, chi_p = y
        
        # Stop integration just before z=0 to avoid numerical issues
        if z < 1e-9:
            return [chi_p, 0.0]

        f_z = f(z, lgb)
        fp_z = fp(z, lgb)
        
        # Avoid division by zero at z=1, though we start integration slightly away
        if abs(f_z) < 1e-12:
            return [0.0, 0.0]

        # EOM term 1: Involves the first derivative χ'
        term1 = -(fp_z / f_z - 1.0 / z) * chi_p
        
        # EOM term 2: Involves the function χ itself.
        # This includes the effect of the chemical potential μ.
        # The term phi_sq_term/f_z is what drives the condensation.
        phi_sq_term = (q_sq_val * mu**2 * (1 - z**2)**2)
        term2_coeff = (m_sq_val - phi_sq_term / f_z) / (z**2 * f_z)
        term2 = term2_coeff * chi

        d_chi_p_dz = term1 + term2
        
        return [chi_p, d_chi_p_dz]

    # --- Define the shooting function ---
    def shoot(mu):
        """
        Solves the ODE for a given μ and returns the value of χ at the boundary.
        The root of this function is the critical chemical potential μ_c.
        """
        # Integrate from near the horizon (z=1) to near the boundary (z=0)
        z_start = 1.0 - 1e-6
        z_end = 1e-5
        
        # Set initial conditions at the horizon for a regular solution
        # χ(1)=1 (by normalization), χ'(1)=0
        y0 = [1.0, 0.0]
        
        # Numerically solve the initial value problem
        sol = solve_ivp(
            fun=ode_system, 
            t_span=[z_start, z_end], 
            y0=y0, 
            args=(mu, lambda_gb, m_sq, q_sq)
        )
        
        # The boundary condition at z=0 is χ=0. We return the final integrated value.
        final_chi = sol.y[0][-1]
        return final_chi

    # --- Find the root to get the critical potential ---
    # We search for the root in a physically reasonable interval, e.g., [1.0, 2.0].
    # By testing, we find the root is between 1.2 and 1.4.
    result = root_scalar(shoot, bracket=[1.2, 1.4], method='brentq')
    mu_c = result.root

    print("The final calculation determines the critical chemical potential, μ_c.")
    print(f"The input parameters used were:")
    print(f"  - Gauss-Bonnet coupling (λ_GB): {lambda_gb}")
    print(f"  - Scalar field mass squared (m²): {m_sq}")
    print(f"  - Scalar field charge squared (q²): {q_sq}")
    print("\nThe computed value is:")
    print(mu_c)

# Execute the calculation
calculate_critical_potential()