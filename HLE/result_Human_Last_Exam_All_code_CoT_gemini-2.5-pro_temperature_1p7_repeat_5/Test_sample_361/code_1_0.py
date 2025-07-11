import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from numpy.polynomial.polynomial import polyfit, polyval

def solve_critical_mu():
    """
    Calculates the critical chemical potential for scalar condensation in a 5D
    Einstein-Gauss-Bonnet holographic model.
    """
    # --- Model Parameters ---
    lgb = 0.1  # Gauss-Bonnet coupling
    m2 = -3.0  # Scalar field mass squared (m^2 * L^2 = -3)
    q = 1.0    # Scalar field charge
    d = 5      # Spacetime dimension

    # In d=5 (d-1=4 dim boundary), the radial coordinate is rho,
    # but z is more common. EOMs are different depending on coordinates.
    # The Maxwell EOM for A_t(z) is phi'' + (g_'/g)phi'=0
    # In standard z coordinates (AdS5), it is phi'' - (1/z)phi' = 0
    # so phi(z) = mu - rho * z^2. Let's use this form.
    # This comes from Maxwell eq: d/dz( (L/z) * phi' ) = 0.
    
    # Background metric function f(z)
    def f(z, zh):
        if z >= zh:
            return 0.0
        # This formula assumes L=1
        radicand = 1.0 - 4.0 * lgb * (1.0 - (z / zh)**4)
        if radicand < 0:
            # Should not happen for z < zh and lgb <= 1/4
            return 1.0 
        return (1.0 - np.sqrt(radicand)) / (2.0 * lgb)

    # Derivative of f(z)
    def f_prime(z, zh):
        if z >= zh:
             # Regularized value at horizon z=zh is -4/zh, but can be singular for numerics
             return -4.0/zh
        radicand = 1.0 - 4.0 * lgb * (1.0 - (z / zh)**4)
        if radicand <= 0:
            return -4.0/zh
        return -4.0 * z**3 / (zh**4 * np.sqrt(radicand))

    # This defines the ODE system for the scalar field Psi(z)
    # Y = [Psi, Psi']
    def odesystem(z, Y, mu, zh):
        Psi, dPsi = Y
        
        # We need values of f(z) and its derivative
        f_val = f(z, zh)
        if f_val == 0.0: # Avoid division by zero at the horizon
            return [0, 0]
        fp_val = f_prime(z, zh)
        
        # Gauge field potential phi(z)
        phi_val = mu * (1.0 - (z / zh)**2)
        
        # Equation for Psi:
        # Psi'' + (f'/f - 3/z)Psi' + (q^2*phi^2/f^2 - m^2/(z^2*f))Psi = 0
        d2Psi = -(fp_val / f_val - (d-2)/z) * dPsi - (q**2 * phi_val**2 / f_val**2 - m2 / (z**2 * f_val)) * Psi
        return [dPsi, d2Psi]

    # This function shoots from the horizon to the boundary and extracts the source term
    def shoot(mu, zh, z_bdy=1e-5):
        # Start integration very close to the horizon
        z_start = zh * (1.0 - 1e-7)
        
        # Set initial conditions for Psi at the horizon
        # Regularity requires Psi'(zh) = m^2/(zh^2*f'(zh)) * Psi(zh)
        f_prime_h = -4.0 / zh 
        psi_h = 1.0
        dpsi_h = (m2 / (zh * f_prime_h)) * psi_h

        # Solve the IVP
        sol = solve_ivp(odesystem, [z_start, z_bdy], [psi_h, dpsi_h], 
                        args=(mu, zh), dense_output=True, atol=1e-9, rtol=1e-9)
        
        # Extract solution at the boundary
        psi_bdy, dpsi_bdy = sol.sol(z_bdy)

        # Asymptotic behavior at z->0 is Psi(z) = C1*z^(d-3-nu) + C2*z^nu
        # Here d=5, nu=sqrt((d-2)^2+m^2)=sqrt(4+m^2)=1
        # Psi(z) ~ psi_1 * z + psi_2 * z^3
        # psi_1 is the source, psi_2 is the condensate. We want psi_1 = 0.
        # From psi(z) and psi'(z), we can find psi_1 and psi_2
        # psi_bdy = psi_1 * z_bdy + psi_2 * z_bdy**3
        # dpsi_bdy = psi_1 + 3 * psi_2 * z_bdy**2
        psi_1 = (3 * psi_bdy - z_bdy * dpsi_bdy) / (2 * z_bdy)
        
        return psi_1
    
    # --- Main Calculation ---
    zh_values = np.array([10.0, 15.0, 20.0, 25.0])
    mu_values = []
    
    print("Calculating critical chemical potential for different horizon radii (zh):")
    # A wide bracket for the root-finder. The value should be of order 1.
    mu_bracket = [1.0, 15.0]
    
    for zh in zh_values:
        # Find mu that makes the source term vanish
        # The lambda function is what we want to find the root of.
        try:
            mu_c_zh = brentq(lambda mu: shoot(mu, zh), mu_bracket[0], mu_bracket[1], xtol=1e-6)
            mu_values.append(mu_c_zh)
            print(f"  For zh = {zh:>4.1f} (T ~ {1/(np.pi*zh):.4f}), found mu_c = {mu_c_zh:.8f}")
            # Update bracket for next iteration to speed up search
            mu_bracket[0] = mu_c_zh * 0.8
            mu_bracket[1] = mu_c_zh * 1.2
        except ValueError:
            print(f"  Root finding failed for zh = {zh}. The bracket might be wrong.")
            # If it fails, we cannot proceed for this point.
            pass
            
    # Extrapolate to zh -> infinity (T -> 0)
    # We fit mu(zh) to a polynomial in (1/zh) and find the intercept.
    # A linear fit is often sufficient for large zh.
    inv_zh = 1.0 / np.array(zh_values)
    # Fit mu = P(1/zh) = c0 + c1*(1/zh)
    coeffs = polyfit(inv_zh, mu_values, 1)
    mu_c_T0 = coeffs[0]
    
    print("\n--- Extrapolation to Zero Temperature ---")
    print(f"Fit function: mu(1/zh) = c0 + c1 * (1/zh)")
    print(f"Fit results for mu at different horizon radii ({zh_values}):")
    # Reconstructing the equation with the found values
    final_equation_parts = []
    for i, (zh, mu) in enumerate(zip(zh_values, mu_values)):
        part = f"{mu:.4f} = {coeffs[0]:.4f} + {coeffs[1]:.4f} * (1/{zh:.1f})"
        final_equation_parts.append(part)
        print(f"  Data point {i+1}: {part}")
        
    print("\nThe critical chemical potential at T=0 is the intercept c0.")
    print(f"μ_c(T=0) = {mu_c_T0:.8f}")

    return mu_c_T0

# Run the calculation and get the final answer
critical_mu = solve_critical_mu()

# The final response should be just the value
print(f"\nFinal answer obtained from extrapolation: {critical_mu}")

<<<8.78310378>>>