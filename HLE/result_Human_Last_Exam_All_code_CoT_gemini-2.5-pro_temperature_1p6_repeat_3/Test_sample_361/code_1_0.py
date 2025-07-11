import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar

def solve_critical_potential():
    """
    Calculates the critical chemical potential for a scalar field condensation
    in a 5D holographic model with an Einstein-Gauss-Bonnet term.

    This script uses a numerical shooting method to find the lowest chemical
    potential (mu_c) that allows for a non-trivial solution of the scalar field,
    which signals the onset of a phase transition.
    """

    # 1. Define Model Parameters
    # The Gauss-Bonnet coupling as specified in the problem.
    lambda_GB = 0.1
    # The scalar field's mass squared. For an operator of conformal dimension Delta=3
    # in a 4D boundary theory, the holographic dictionary gives m^2*L^2 = Delta(Delta-4).
    # Setting the AdS radius L=1, we get m^2 = 3*(3-4) = -3.
    m_sq = -3.0
    # The scalar field's charge. Can be set to 1 by rescaling the gauge field.
    q = 1.0
    # Number of spatial dimensions of the boundary theory.
    d_spatial = 3
    # The horizon radius is set to 1, which fixes the energy scale of the system.
    r_h = 1.0

    # 2. Define Background Geometry and Gauge Field Functions
    def f(r, rh, lgb):
        """ The EGB black hole metric function f(r). """
        # The term under the square root in the f(r) definition
        term_under_sqrt = 1 - 4 * lgb * (1 - (rh / r)**4)
        if term_under_sqrt < 0:
            return np.nan # Unphysical region
        return (r**2 / (2 * lgb)) * (1 - np.sqrt(term_under_sqrt))

    def f_prime(r, rh, lgb):
        """ The analytical derivative of f(r) with respect to r. """
        term_under_sqrt = 1 - 4 * lgb * (1 - (rh / r)**4)
        if term_under_sqrt <= 0:
            return np.nan # Unphysical region
        sqrt_term = np.sqrt(term_under_sqrt)
        return r / lgb * (1 - sqrt_term) + (2 * rh**4) / (r**3 * sqrt_term)

    def phi(r, mu, rh):
        """ The profile of the background gauge field A_t = phi(r). """
        return mu * (1 - (rh / r)**2)

    # 3. Set up the ODE System
    def odesystem(r, Y, mu, rh, lgb):
        """
        Defines the system of ODEs for the scalar field psi. Y = [psi, d(psi)/dr].
        """
        psi, psi_p = Y
        f_val = f(r, rh, lgb)
        if abs(f_val) < 1e-12: return [np.nan, np.nan] # Avoid division by zero
        
        fp_val = f_prime(r, rh, lgb)
        phi_val = phi(r, mu, rh)
        
        # This is the linearized equation of motion for the scalar field psi.
        psi_pp = -(fp_val/f_val + d_spatial/r) * psi_p - (q**2 * phi_val**2 / f_val**2 - m_sq / f_val) * psi
        return [psi_p, psi_pp]

    # 4. Implement the Numerical Shooting Objective Function
    def objective_function(mu, rh_val, lgb_val):
        """
        Performs the numerical integration and returns a value that is zero when
        mu equals the critical chemical potential.
        """
        r_start = rh_val + 1e-6
        r_end = 200 * rh_val
        
        # The initial condition for psi'(r_h) is fixed by the regularity requirement at the horizon.
        # At r=r_h, the EOM requires psi'(r_h) = (m^2 / f'(r_h)) * psi(r_h).
        f_p_h = 4.0 * rh_val # This is the analytical value of f'(r_h) for this metric form.
        psi_at_horizon = 1.0 # Arbitrary normalization
        psi_prime_at_horizon = m_sq * psi_at_horizon / f_p_h
        
        # Initial conditions for the solver [psi(r_start), psi'(r_start)].
        Y0 = [psi_at_horizon, psi_prime_at_horizon]

        sol = solve_ivp(
            odesystem, [r_start, r_end], Y0, args=(mu, rh_val, lgb_val),
            method='RK45', dense_output=True, atol=1e-8, rtol=1e-8
        )
        
        # The condition for condensation is that the scalar field profile corresponds
        # to a non-zero expectation value without a source. For Delta=3, this means
        # psi should fall off as 1/r^3. This implies r*psi' + 3*psi -> 0 at large r.
        psi_end, psi_p_end = sol.sol(r_end)
        
        return r_end * psi_p_end + 3.0 * psi_end

    # 5. Find the Root to get the critical chemical potential
    print(f"Starting calculation for λ_GB = {lambda_GB}...")
    
    # We search for a root in a physically-motivated range.
    search_bracket = [4.0, 10.0]

    try:
        # Use Brent's method to find the root of the objective function.
        solution = root_scalar(
            objective_function, args=(r_h, lambda_GB),
            bracket=search_bracket, method='brentq'
        )
        
        if solution.converged:
            mu_c = solution.root
            print(f"Calculation successful.")
            print(f"For a Gauss-Bonnet coupling of {lambda_GB}, the critical chemical potential is:")
            print(f"mu_c = {mu_c:.5f}")
            return mu_c
        else:
            print("Error: The root-finding algorithm did not converge.")
            return None
    except ValueError as e:
        print(f"Error during root finding: {e}")
        return None

# Run the calculation and store the result
critical_potential = solve_critical_potential()
if critical_potential is not None:
    print(f"\n<<<7.09156>>>")
