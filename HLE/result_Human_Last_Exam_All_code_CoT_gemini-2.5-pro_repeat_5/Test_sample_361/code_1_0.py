import numpy as np
from scipy.integrate import solve_ivp
import warnings

def solve_critical_potential():
    """
    Calculates the critical chemical potential for scalar condensation in a D3/D7 holographic model
    with an Einstein-Gauss-Bonnet background.
    """
    # Suppress warnings that can occur during the numerical integration of stiff ODEs.
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    # 1. Define physical constants of the model
    LAMBDA_GB = 0.1  # Gauss-Bonnet coupling
    RH = 1.0         # Horizon radius, which sets the temperature scale (T = rh/pi)
    M2 = -3.0        # Scalar field mass squared in AdS units (m^2 * L^2)

    # 2. Define the EGB blackening factor f(r) and its derivative f'(r)
    # This is a numerically stable form of f(r) to avoid issues near the horizon.
    def f(r, lambda_gb):
        x = 1.0 - (RH**4 / r**4)
        if abs(x) < 1e-14:
            return 0.0
        # This expression is equivalent to (r**2 / (2*lambda_gb)) * (1 - sqrt(1 - 4*lambda_gb*x))
        return (2.0 * r**2 * x) / (1.0 + np.sqrt(1.0 - 4.0 * lambda_gb * x))

    def dfdr(r, lambda_gb):
        # The derivative at the horizon has a well-defined limit: f'(rh) = 4*rh
        if abs(r - RH) < 1e-7:
            return 4.0 * RH

        g = 1.0 - 4.0 * lambda_gb * (1.0 - (RH**4 / r**4))
        if g < 1e-12:
            return 4.0 * RH # Fallback to the limit value

        sqrt_g = np.sqrt(g)
        g_prime = -4.0 * lambda_gb * (4.0 * RH**4 / r**5)
        
        # Standard derivative calculation
        term1 = (r / lambda_gb) * (1.0 - sqrt_g)
        term2 = (r**2 / (4.0 * lambda_gb * sqrt_g)) * g_prime
        return term1 - term2

    # 3. Define the system of ODEs for the scalar field perturbation Phi
    def ode_system(r, y, mu, lambda_gb):
        phi, dphi_dr = y

        fr = f(r, lambda_gb)
        if abs(fr) < 1e-12:
            return [dphi_dr, 0.0]
        
        dfr = dfdr(r, lambda_gb)

        # Background gauge field A_t(r) corresponds to the chemical potential
        A_t = mu * (1.0 - (RH**2 / r**2))

        # The scalar field equation: Phi'' + P(r)*Phi' + Q(r)*Phi = 0
        P_r = (3.0 / r) + (dfr / fr)
        Q_r = (A_t**2 / fr**2) + (M2 / fr)
        
        ddphi_dr2 = -P_r * dphi_dr - Q_r * phi
        return [dphi_dr, ddphi_dr2]

    # 4. Shooting function: solves the ODE for a given mu and returns Phi at r_end
    def get_phi_at_infinity(mu, r_start, r_end):
        # Initial conditions from regularity at the horizon: Phi(rh)=1, Phi'(rh)=-m^2/f'(rh)
        y0 = [1.0, -M2 / (4.0 * RH)]
        
        sol = solve_ivp(
            ode_system, [r_start, r_end], y0,
            args=(mu, LAMBDA_GB), method='RK45', dense_output=True,
            atol=1e-8, rtol=1e-8
        )
        
        if sol.status != 0:
            return np.inf # Return infinity to indicate integration failure
        return sol.sol(r_end)[0]

    # 5. Bisection search to find mu_c where get_phi_at_infinity(mu_c) = 0
    r_start = RH + 1e-6
    r_end = 300.0

    # Define a bracket [mu_low, mu_high] that should contain the root
    mu_low, mu_high = 0.5, 1.5

    val_low = get_phi_at_infinity(mu_low, r_start, r_end)
    val_high = get_phi_at_infinity(mu_high, r_start, r_end)

    if np.sign(val_low) == np.sign(val_high):
        print("Error: The bisection method bracket is invalid. The function does not change sign.")
        return

    # Perform bisection
    for _ in range(50):
        mu_mid = (mu_low + mu_high) / 2.0
        if mu_mid == mu_low or mu_mid == mu_high:
            break # Precision limit reached
        
        val_mid = get_phi_at_infinity(mu_mid, r_start, r_end)

        if np.sign(val_mid) == np.sign(val_low):
            mu_low = mu_mid
        else:
            mu_high = mu_mid
    
    mu_c = (mu_low + mu_high) / 2.0

    # Print the final result in a descriptive way
    print("In the context of the specified D3/D7 holographic model:")
    print(f"Gauss-Bonnet coupling (lambda_GB): {LAMBDA_GB}")
    print(f"Scalar field mass squared (m^2): {M2}")
    print(f"Horizon radius (r_h): {RH} (This sets the temperature T = r_h/pi)")
    print("\nThe equation for the scalar perturbation Phi(r) takes the form:")
    print("Phi''(r) + P(r)*Phi'(r) + Q(r)*Phi(r) = 0, where")
    print("Q(r) = (mu^2 * (1 - r_h^2/r^2)^2 / f(r)^2) + m^2/f(r)")
    print("\nThe critical chemical potential 'mu_c' is found when this equation has a normalizable solution.")
    print("\nResult:")
    print(f"The critical chemical potential is mu_c = {mu_c:.6f}")
    
    return mu_c

# Execute the calculation
critical_mu = solve_critical_potential()
# Final answer in the required format
if critical_mu is not None:
    print(f"\n<<<{critical_mu:.6f}>>>")
