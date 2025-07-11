import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def solve_critical_potential():
    """
    Calculates the critical chemical potential for a holographic superconductor
    in 5D Einstein-Gauss-Bonnet gravity.
    """
    # 1. Define the physical parameters of the model
    # Gauss-Bonnet coupling constant
    LAMBDA_GB = 0.1
    # Mass squared of the scalar field. For the operator of dimension 3 in the D3/D7 model, m^2 L^2 = -3.
    M2 = -3.0
    # Charge of the scalar field, can be set to 1 without loss of generality.
    Q = 1.0

    # 2. Calculate the asymptotic properties of the solution
    # The EGB background behaves like a different AdS space at infinity,
    # with an effective radius determined by k^2.
    k2 = (1 - np.sqrt(1 - 4 * LAMBDA_GB)) / (2 * LAMBDA_GB)
    
    # The scalar field at infinity behaves as psi ~ c1*r^nu_plus + c2*r^nu_minus.
    # nu_plus corresponds to the source, nu_minus to the expectation value.
    # We need to find mu_c for which the source c1 is zero.
    falloff_sqrt = np.sqrt(4 + M2 / k2)
    nu_plus = -2 + falloff_sqrt  # Slower falloff (source)
    nu_minus = -2 - falloff_sqrt # Faster falloff (VEV)

    # 3. Define the metric functions for the EGB black hole
    def f(r):
        """Blackening function f(r) for the EGB metric."""
        # Add a small epsilon for numerical stability near r=1
        if np.isclose(r, 1.0):
            return 0.0
        arg_sqrt = 1 - 4 * LAMBDA_GB * (1 - 1 / (r**4 + 1e-15))
        # This argument should not be negative for the allowed range of lambda_gb
        if arg_sqrt < 0:
            arg_sqrt = 0
        return (r**2 / (2 * LAMBDA_GB)) * (1 - np.sqrt(arg_sqrt))

    def f_prime(r):
        """Derivative of the blackening function, f'(r)."""
        if np.isclose(r, 1.0):
            return 4.0 # Analytic value at r_h=1
        arg_sqrt = 1 - 4 * LAMBDA_GB * (1 - 1 / r**4)
        if arg_sqrt <= 0:
            return np.inf
        sqrt_term = np.sqrt(arg_sqrt)
        term1 = r / LAMBDA_GB * (1 - sqrt_term)
        term2 = (4 * r**-3) / sqrt_term
        return term1 + term2

    # 4. Define the system of ODEs for the scalar field psi
    def ode_system(r, y, mu):
        """
        Converts the second-order ODE for psi into a system of two first-order ODEs.
        y = [psi, dpsi/dr]
        """
        psi, dpsi = y
        
        # In the normal phase, the gauge field has a simple analytic solution
        phi = mu * (1 - 1/r**2)
        fr = f(r)
        
        # Avoid division by zero precisely at the horizon
        if np.isclose(fr, 0.0):
            return [0, 0] # Should not be evaluated here by the solver

        # This is the main equation: d2psi/dr2 = - (f'/f + 3/r) dpsi + (m^2 - q^2*phi^2/f)/f * psi
        d2psi = -(f_prime(r)/fr + 3/r) * dpsi + (M2 - Q**2 * phi**2 / fr) / fr * psi
        
        return [dpsi, d2psi]

    # 5. Define the objective function for the root-finder
    def objective_function(mu):
        """
        This function returns a value that is zero when the boundary condition at infinity is met.
        We will find the root of this function to get mu_c.
        """
        # Integration range
        r_start = 1.0 + 1e-6
        r_end = 200.0

        # Initial conditions at the horizon (r=1) that ensure regularity
        f_prime_at_1 = 4.0 
        psi_0 = 1.0 # The overall scale is arbitrary for this linear ODE
        dpsi_0 = M2 / f_prime_at_1 * psi_0

        # Solve the ODE system using a Runge-Kutta method
        sol = solve_ivp(
            fun=lambda r, y: ode_system(r, y, mu),
            t_span=[r_start, r_end],
            y0=[psi_0, dpsi_0],
            method='RK45',
            dense_output=True,
            atol=1e-9, rtol=1e-9
        )
        
        # If the solver fails, return a large value to guide the root finder away
        if sol.status != 0:
            return np.sign(mu - 5.0) * 1e10

        # Extract the solution at the endpoint r_end
        psi_end, dpsi_end = sol.y[:, -1]

        # The target for the root finder. It is zero if psi ~ r^nu_minus at infinity.
        # This condition is equivalent to setting the source coefficient to zero.
        target = dpsi_end * r_end - nu_minus * psi_end
        return target

    # 6. Find the root of the objective function to get the critical potential mu_c
    try:
        # We search for the root in a physically plausible interval, e.g., [3.0, 8.0]
        mu_c = brentq(objective_function, a=3.0, b=8.0, xtol=1e-6)
        
        # 7. Print the final results, including the parameters used
        print("Calculation of the critical chemical potential for scalar condensation.")
        print("-" * 70)
        print(f"Model Parameters:")
        print(f"  Gauss-Bonnet coupling (lambda_GB) = {LAMBDA_GB}")
        print(f"  Scalar field mass squared (m^2)   = {M2}")
        print(f"  Scalar field charge (q)           = {Q}")
        print("-" * 70)
        print("Result:")
        print(f"The critical chemical potential is mu_c = {mu_c:.4f}")
        
    except ValueError:
        print("Error: Could not find the critical potential in the given range [3.0, 8.0].")
        print("The root may not be bracketed by these values.")

if __name__ == '__main__':
    solve_critical_potential()