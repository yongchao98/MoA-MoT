import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

# Define the model parameter
LAMBDA_GB = 0.1

# Define the metric function f(z) and its derivative f'(z)
# z_h is set to 1.
def f(z):
    """EGB black hole metric function f(z)."""
    if z >= 1:
        return 0.0
    val_in_sqrt = 1 - 4 * LAMBDA_GB * (1 - z**4)
    if val_in_sqrt < 0:
        return np.nan # Should not happen for z in [0,1] and LAMBDA_GB < 1/4
    return (1 - np.sqrt(val_in_sqrt)) / (2 * LAMBDA_GB)

def f_prime(z):
    """Derivative of the metric function f(z) with respect to z."""
    if z == 0:
        return 0.0
    if z >= 1:
        # A more stable evaluation near z=1, even though it's formally singular
        return -4 / np.sqrt(1 - 4 * LAMBDA_GB)
        
    val_in_sqrt = 1 - 4 * LAMBDA_GB * (1 - z**4)
    if val_in_sqrt <= 0:
         return np.nan
    return (8 * LAMBDA_GB * z**3) / (2 * LAMBDA_GB * np.sqrt(val_in_sqrt))

def ode_system(z, y, mu):
    """
    Defines the ODE system for the scalar field psi.
    y[0] = psi, y[1] = psi'
    """
    psi, psi_prime = y
    
    # Avoid singularities at z=0 and z=1 by returning 0 for derivatives
    # The solver will handle small steps near the endpoints.
    if z == 0 or z == 1:
        return [0, 0]
        
    f_val = f(z)
    fp_val = f_prime(z)
    
    # The ODE: psi'' + P(z)psi' + Q(z)psi = 0
    # psi'' = -P(z)psi' - Q(z)psi
    # P(z) = f'/f - 1/z
    # Q(z) = mu^2 * (1-z^2)^2 / f^2
    
    d_psi_dz = psi_prime
    d_psi_prime_dz = -(fp_val / f_val - 1/z) * psi_prime - (mu**2 * (1 - z**2)**2 / f_val**2) * psi
    
    return [d_psi_dz, d_psi_prime_dz]

def shoot(mu):
    """
    Shooting function. Solves the IVP for a given mu and returns the
    residual at the boundary z=0. The root of this function is mu_c.
    """
    # Small number to avoid singularities at z=1 and z=0
    epsilon = 1e-6
    z_start = 1.0 - epsilon
    z_end = epsilon

    # Calculate f'(1)
    f1_val = np.abs(f_prime(1.0))

    # Boundary conditions near the horizon (z=1) from Frobenius method
    # We set psi(1) = 1 (normalization)
    # This implies psi'(1) = -4*mu^2/f1^2
    psi_prime_at_1 = -4 * mu**2 / f1_val**2
    
    # Initial conditions at z = 1 - epsilon
    # psi(1-eps) ~ psi(1) + psi'(1)*(z-1) = 1 + psi'(1)*(-eps)
    psi_start = 1.0 - psi_prime_at_1 * epsilon
    psi_prime_start = psi_prime_at_1
    y0 = [psi_start, psi_prime_start]

    # Solve the ODE from near the horizon to near the boundary
    sol = solve_ivp(
        fun=ode_system,
        t_span=[z_start, z_end],
        y0=y0,
        args=(mu,),
        dense_output=True,
        method='RK45'
    )
    
    # Extract the solution at the endpoint z=epsilon
    psi_end, psi_prime_end = sol.y[:,-1]
    
    # The condition for a normalizable solution at z=0 is that the
    # coefficient of the non-normalizable mode vanishes. This translates to:
    # psi(epsilon) - epsilon * psi'(epsilon) = 0
    residual = psi_end - epsilon * psi_prime_end
    
    return residual

# --- Main execution ---
# Bracket for the root finding algorithm
mu_min = 1.0
mu_max = 3.0

try:
    # Find the critical chemical potential mu_c using a root finder
    mu_c = brentq(shoot, mu_min, mu_max, xtol=1e-8)
    
    # Output the final result in a descriptive sentence
    final_equation_string = f"For a Gauss-Bonnet coupling of lambda_GB = {LAMBDA_GB}, the equation mu_c = {mu_c:.8f} gives the critical chemical potential."
    print(final_equation_string)

except ValueError:
    print(f"Failed to find a root in the interval [{mu_min}, {mu_max}].")
    print("The shoot function might not be crossing zero in this interval.")
    print(f"Value at mu_min={mu_min}: {shoot(mu_min)}")
    print(f"Value at mu_max={mu_max}: {shoot(mu_max)}")
