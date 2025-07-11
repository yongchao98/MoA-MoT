import numpy as np
from scipy.integrate import dblquad

def calculate_hopf_charge():
    """
    Calculates the Hopf charge for the given vector field using a 
    simplified 2D integral derived from the Whitehead formula.
    """
    PI = np.pi

    # Define the function G(rho, z)
    def G(rho, z):
        # To avoid division by zero or sqrt of negative numbers, handle small values
        rho = max(rho, 1e-9)
        r2_sq = (rho**2 - 0.5)**2 + z**2
        if r2_sq < 1e-16:
            # This is the center of the skyrmion ring, where G = pi
            return PI
        r2 = np.sqrt(r2_sq)
        return PI * np.exp(-10.0 * r2)

    # Define the partial derivative of G with respect to rho
    def dG_drho(rho, z):
        # By symmetry, the derivative is 0 on the z-axis (rho=0)
        if rho < 1e-9:
            return 0.0
        
        r2_sq = (rho**2 - 0.5)**2 + z**2
        # Avoid division by zero at the singularity
        if r2_sq < 1e-16:
            # The integrand is handled separately at the singularity, 
            # but the derivative itself is undefined here.
            # We return NaN and the integrand function will handle it.
            return np.nan

        r2 = np.sqrt(r2_sq)
        exp_term = np.exp(-10.0 * r2)
        
        # Chain rule: dG/drho = dG/dr2 * dr2/drho
        dG_dr2 = -10.0 * PI * exp_term
        dr2_drho = (2.0 * rho * (rho**2 - 0.5)) / r2
        
        return dG_dr2 * dr2_drho

    # Define the integrand for the double integral
    # Note: dblquad expects the function signature func(y, x),
    # so we map y -> z and x -> rho.
    def integrand(z, rho):
        # The integrand is singular at the point (rho=sqrt(0.5), z=0),
        # but the singularity is integrable. For numerical stability,
        # we can return 0 right at the singularity, as its contribution
        # to the integral over an infinitesimal area is zero.
        r2_sq = (rho**2 - 0.5)**2 + z**2
        if r2_sq < 1e-16:
            return 0.0
            
        g_val = G(rho, z)
        dG_val = dG_drho(rho, z)
        
        return (1.0 - np.cos(g_val)) * dG_val

    # Perform the numerical integration.
    # The integrand decays very fast, so we can use finite limits.
    rho_limit = 5.0
    z_limit = 5.0
    
    # The dblquad function returns the result and an estimated error.
    integral_value, error = dblquad(integrand, 0, rho_limit, -z_limit, z_limit)

    # Calculate the Hopf charge using the formula H = -1/(2*pi) * integral
    hopf_charge = -1.0 / (2.0 * PI) * integral_value

    print(f"The numerical value of the integral is: {integral_value}")
    print(f"The estimated error of the integral is: {error}")
    print("\nThe Hopf charge H is calculated as:")
    print(f"H = -1/(2*pi) * ({integral_value:.6f})")
    print(f"H = {hopf_charge:.6f}")
    print(f"\nAs a topological invariant, the Hopf charge is an integer.")
    print(f"The calculated value rounded to the nearest integer is: {round(hopf_charge)}")

calculate_hopf_charge()