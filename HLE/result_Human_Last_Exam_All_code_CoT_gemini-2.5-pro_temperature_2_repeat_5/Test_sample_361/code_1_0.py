import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import warnings

# Suppress warnings that may arise during integration near singularities
warnings.filterwarnings("ignore", category=RuntimeWarning)

# --- Model Parameters ---
# Gauss-Bonnet coupling constant
LAMBDA_GB = 0.1
# Squared mass of the scalar field for a Delta=3 operator
M2 = -3.0
# Horizon position, set to 1 for simplicity
ZH = 1.0

def f(z, lambda_gb):
    """
    Defines the blackening function f(z) for the EGB metric.
    """
    # Clamp z to avoid numerical issues at z=1
    z_safe = np.minimum(z, ZH - 1e-9)
    # Argument of the square root
    sqrt_arg = 1.0 - 4.0 * lambda_gb * (1.0 - (z_safe/ZH)**4)
    if sqrt_arg < 0:
        # This signifies a causality violation in the bulk, but our parameters avoid this.
        return np.nan
    return (1.0 / (2.0 * lambda_gb)) * (1.0 - np.sqrt(sqrt_arg))

def f_prime(z, lambda_gb):
    """
    Defines the derivative of the blackening function, f'(z).
    """
    denominator = np.sqrt(1.0 - 4.0 * lambda_gb * (1.0 - (z/ZH)**4))
    if denominator == 0:
        return np.inf
    return 4.0 * z**3 / (ZH**4 * denominator)

def get_boundary_exponent(lambda_gb):
    """
    Calculates the exponent for the condensate operator at the boundary z=0.
    """
    if lambda_gb == 0:
        # Standard AdS case
        return 3.0
    else:
        # Value of f(z) at the boundary z=0
        f0 = (1.0 / (2.0 * lambda_gb)) * (1.0 - np.sqrt(1.0 - 4.0 * lambda_gb))
        # Effective mass at the boundary
        m_eff_sq = M2 / f0
        # The exponent of the operator condensate term (larger root of indicial eq)
        alpha_plus = 2.0 + np.sqrt(4.0 + m_eff_sq)
        return alpha_plus

def ode_system(z, y, mu, lambda_gb):
    """
    Defines the system of first-order ODEs for the scalar field psi.
    y = [psi, dpsi/dz]
    """
    psi, dpsi_dz = y
    
    # Avoid singularity at z=0
    if z < 1e-9:
        return [dpsi_dz, 0]

    f_val = f(z, lambda_gb)
    fp_val = f_prime(z, lambda_gb)
    
    # Coefficients of the ODE: psi'' + P(z)psi' + Q(z)psi = 0
    P_z = (fp_val / f_val) - (3.0 / z)
    Q_z = ((mu**2 * (1.0 - (z/ZH)**2)**2) / f_val**2) - (M2 / (z**2 * f_val))
    
    d2psi_dz2 = -P_z * dpsi_dz - Q_z * psi
    return [dpsi_dz, d2psi_dz2]

def shoot(mu, lambda_gb, alpha_plus):
    """
    Shooting method function. For a given mu, integrates the ODE and
    returns a value that is zero when the boundary conditions are met.
    """
    # Integration range starts near the horizon and goes to near the boundary
    z_start = ZH - 1e-6
    z_end = 1e-5
    
    # Initial conditions from regularity at the horizon psi'(z_h) = (m^2/f'(z_h)) * psi(z_h)
    psi_start = 1.0  # Normalization
    psi_prime_start = (M2 / f_prime(ZH, lambda_gb)) * psi_start
    y_start = [psi_start, psi_prime_start]
    
    # Solve the ODE system
    sol = solve_ivp(fun=ode_system, t_span=[z_start, z_end], y0=y_start,
                    args=(mu, lambda_gb), dense_output=True, method='RK45')
    
    # Extract the solution at the endpoint near the boundary
    psi_end, dpsi_dz_end = sol.sol(z_end)
    
    # The condition for the condensate solution (c_src = 0)
    # is that z * psi' = alpha_plus * psi at z -> 0.
    # We want to find mu such that the following expression is zero.
    target_value = z_end * dpsi_dz_end - alpha_plus * psi_end
    return target_value

def find_critical_potential(lambda_gb):
    """
    Finds the critical chemical potential mu_c using a root-finding algorithm.
    """
    alpha_plus = get_boundary_exponent(lambda_gb)

    try:
        # Bracket the root. Based on known physics, the result for lambda_gb=0.1
        # should be in this range.
        mu_bracket = [1.5, 2.0]
        mu_c = brentq(shoot, mu_bracket[0], mu_bracket[1], args=(lambda_gb, alpha_plus), xtol=1e-6)
        return mu_c
    except ValueError:
        print("Root not found in the initial bracket. The bracket may need adjustment.")
        return None

if __name__ == '__main__':
    critical_mu = find_critical_potential(LAMBDA_GB)

    print("This script calculates the critical chemical potential for scalar condensation")
    print("in a holographic D3/D7 model with Einstein-Gauss-Bonnet gravity in the bulk.")
    print("\n--- Input Parameters ---")
    print(f"Gauss-Bonnet coupling (lambda_GB) = {LAMBDA_GB}")
    print(f"Scalar field mass squared (m^2)   = {M2}")

    if critical_mu is not None:
        print("\n--- Result ---")
        print("The numerical calculation yields the following critical chemical potential:")
        # The question asks to output numbers in a final equation. Since no equation template
        # was given, we present the result clearly.
        final_mu_value = f"{critical_mu:.4f}"
        print(f"mu_c = {final_mu_value}")
    else:
        print("\n--- Result ---")
        print("Could not determine the critical chemical potential.")
>>>1.7486