import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar

def solve_for_critical_potential():
    """
    Calculates the critical chemical potential for scalar condensation in a D3/D7
    holographic model with an Einstein-Gauss-Bonnet background.
    """

    # 1. Define constants and parameters
    LAMBDA_GB = 0.1
    # Mass squared for the scalar field dual to the quark bilinear operator
    M2L2 = -3

    # Define a small number for integration bounds to avoid singularities
    epsilon = 1e-6

    # 2. Define metric function f(z) and its derivative f'(z)
    def f(z, lam):
        """EGB black hole metric function f(z) with L=1, z_h=1."""
        if np.isclose(z, 1.0):
            return 0.0
        # Argument of the square root
        sqrt_arg = 1 - 4 * lam * (1 - z**4)
        if sqrt_arg < 0:
            sqrt_arg = 0 # Clamp due to potential floating point inaccuracies
        return (1 - np.sqrt(sqrt_arg)) / (2 * lam)

    def f_prime(z, lam):
        """Derivative of f(z) with respect to z."""
        if np.isclose(z, 1.0):
           # The limit of f'(z) as z->1 is -4, independent of lambda
           return -4.0
        sqrt_arg = 1 - 4 * lam * (1 - z**4)
        if sqrt_arg <= 0:
            sqrt_arg = epsilon # Avoid division by zero
        return (4 * z**3) / np.sqrt(sqrt_arg)

    # 3. Define the system of first-order ODEs from the scalar field EOM
    # The EOM is Φ'' + P(z)Φ' + Q(z)Φ = 0, which we rewrite as a system
    # for Y = [Φ, Φ'].
    def ode_system(z, Y, mu, lam):
        """
        Defines the system of ODEs for the shooting method.
        Y = [Φ, Φ']
        """
        Phi, dPhi = Y
        
        fz = f(z, lam)
        fpz = f_prime(z, lam)

        # Handle numerical singularities at boundaries
        if np.isclose(z, 0.0): z = epsilon
        if np.isclose(fz, 0.0): fz = epsilon**2

        # Coefficients of the ODE: Φ'' + A(z)Φ' + B(z)Φ = 0
        A = fpz / fz - 3 / z
        B = (M2L2) / (z**2 * fz) - (mu**2 * (1 - z**2)**2) / fz**2
        
        # dY/dz = [Φ', Φ'']
        d2Phi = -A * dPhi - B * Phi
        return [dPhi, d2Phi]

    # 4. Define the shooting function to find the source term
    def get_source_term(mu, lam):
        """
        Shoots from the horizon to the boundary for a given mu and returns the
        source term of the scalar field. The critical potential is found when
        this source is zero.
        """
        if mu <= 0:
            return 1.0  # Source must be non-zero for μ=0

        # Integration range from near horizon (z=1) to near boundary (z=0)
        z_span = [1 - epsilon, epsilon]
        
        # Initial conditions: Regularity at the horizon requires Φ to be constant
        # and Φ' to be zero. We can normalize Φ(1)=1.
        y0 = [1.0, 0.0]

        # Solve the ODE
        sol = solve_ivp(
            ode_system, 
            z_span, 
            y0, 
            args=(mu, lam), 
            dense_output=True,
            method='RK45'
        )
        
        # Extract solution at the boundary z=epsilon
        z_b = epsilon
        Phi_b, dPhi_b = sol.sol(z_b)
        
        # Near the boundary z=0, the solution behaves as Φ(z) ≈ m_q*z + c*z^3.
        # m_q is the source, c is the condensate. We need to find mu where m_q=0.
        # We can extract m_q from the values at z_b=epsilon.
        # 2 * m_q = 3 * Φ(z_b)/z_b - Φ'(z_b)
        source = (3 * Phi_b / z_b - dPhi_b) / 2
        return source

    # 5. Find the root of the source function to get the critical potential
    # Based on known physics, μ_c should be of order 1. We test the bracket
    # [0.5, 1.5] to ensure the root is captured.
    try:
        sol = root_scalar(
            lambda mu: get_source_term(mu, LAMBDA_GB), 
            bracket=[0.5, 1.5], 
            method='brentq'
        )
        critical_mu = sol.root
        
        # The prompt requires outputting each number in the final equation.
        # Let's present the result clearly.
        print("The value of the Gauss-Bonnet coupling lambda_GB is:", LAMBDA_GB)
        print("The value of the critical chemical potential mu_c is:", critical_mu)

    except ValueError:
        print("Could not find the root in the specified bracket.")
        print("The source function may not change sign in the given range.")

if __name__ == '__main__':
    solve_for_critical_potential()
