import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def solve_critical_potential():
    """
    Calculates the critical chemical potential for a holographic superconductor
    in 5D Gauss-Bonnet gravity using a shooting method.
    """
    # Problem Parameters
    LAMBDA_GB = 0.1  # Gauss-Bonnet coupling
    M2 = -3.0       # Mass squared of the scalar field
    ZH = 1.0        # Horizon position (sets the energy scale)

    # --- Equation of Motion Details ---
    # The ODE for the scalar field psi(z) is:
    # psi'' + P(z)psi' + Q(z)psi = 0
    # where:
    # P(z) = f'(z)/f(z) - 3/z
    # Q(z) = (mu^2 * (1 - z^2)^2) / f(z)^2 - m^2 / (z^2 * f(z))
    # f(z) = (1/(2*LAMBDA_GB)) * (1 - sqrt(1 - 4*LAMBDA_GB*(1 - (z/ZH)^4)))

    # Define the metric function f(z) and its derivative f'(z)
    def f(z):
        # Prevent floating point issues near the horizon
        if np.isclose(z, ZH):
            return 0.0
        val_in_sqrt = 1.0 - 4.0 * LAMBDA_GB * (1.0 - (z/ZH)**4)
        return (1.0 / (2.0 * LAMBDA_GB)) * (1.0 - np.sqrt(val_in_sqrt))

    def f_prime(z):
        if z == 0.0:
            return 0.0
        denominator = np.sqrt(1.0 - 4.0 * LAMBDA_GB * (1.0 - (z/ZH)**4))
        # f'(z_h) = -4/z_h is finite
        return -(4.0 * z**3 / ZH**4) / denominator

    # Define the system of first-order ODEs for the solver
    def ode_system(z, y, mu):
        psi, chi = y  # y = [psi, psi_prime]
        
        fz = f(z)
        # Avoid division by zero at the horizon
        if fz == 0:
            return [chi, 0] 

        fz_p = f_prime(z)
        
        P_z = fz_p / fz - 3.0 / z
        Q_z = (mu**2 * (1 - (z/ZH)**2)**2) / fz**2 - M2 / (z**2 * fz)
        
        psi_double_prime = -P_z * chi - Q_z * psi
        return [chi, psi_double_prime]

    # Power law exponents at the boundary z->0 (for extracting the source term)
    f0 = f(1e-9)
    nu_plus = 2 + np.sqrt(4 + M2/f0)
    
    # Shooting function: integrates the ODE for a given mu and returns a value
    # that is zero when the boundary condition at z=0 is met.
    def shoot(mu):
        z_start = ZH - 1e-5
        z_end = 1e-6

        # Initial conditions near the horizon from analytic expansion
        # psi(z) ~ (ZH - z)^r_h where r_h = sqrt(-m^2)/2
        r_h = np.sqrt(-M2) / 2.0
        epsilon = ZH - z_start
        psi_init = epsilon**r_h
        chi_init = -r_h * epsilon**(r_h - 1)
        y_init = [psi_init, chi_init]

        # Solve the ODE
        sol = solve_ivp(
            ode_system, [z_start, z_end], y_init, args=(mu,),
            dense_output=True, method='RK45', atol=1e-8, rtol=1e-8
        )

        psi_end = sol.sol(z_end)[0]
        chi_end = sol.sol(z_end)[1]

        # Target function: c1 ~ z*psi' - nu+ * psi
        # This function is zero when the coefficient of the source term (c1) vanishes.
        target = z_end * chi_end - nu_plus * psi_end
        return target

    # Find the root of the shooting function to get the critical chemical potential
    # Based on literature, the value for mu_c is expected to be in the range [1, 2] for z_h=1.
    mu_c = brentq(shoot, 1.0, 2.0)
    
    # Temperature of the black hole
    T = abs(f_prime(ZH)) / (4 * np.pi)

    # Output the final results of the calculation
    print("Calculation of the critical chemical potential (mu_c)")
    print("-------------------------------------------------------")
    print(f"Model parameters:")
    print(f"  Gauss-Bonnet coupling (lambda_GB): {LAMBDA_GB}")
    print(f"  Scalar field mass squared (m^2):   {M2}")
    print("\nThe final equation for the scalar field psi(z) is:")
    print("psi''(z) + (f'(z)/f(z) - 3/z)psi'(z) + (mu_c^2*(1 - z^2)^2/f(z)^2 + 3/(z^2*f(z)))psi(z) = 0")
    print("where f(z) = 5*(1 - sqrt(0.6 + 0.4*z^4))\n")
    print(f"Calculated Results (for z_h=1):")
    print(f"  Critical Chemical Potential (mu_c): {mu_c:.6f}")
    print(f"  Black Hole Temperature (T):         {T:.6f}")
    print(f"  Dimensionless Ratio (mu_c / T):   {mu_c/T:.6f}")

solve_critical_potential()