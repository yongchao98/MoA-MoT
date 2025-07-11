import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def calculate_critical_potential():
    """
    Calculates the critical chemical potential for a holographic superconductor
    in a 5D Einstein-Gauss-Bonnet background using a shooting method.
    """

    # --- 1. Define constants and model parameters ---
    LAMBDA_GB = 0.1
    # For a scalar dual to an operator of dimension Delta=2 in this setup, m^2*L^2 = -2.
    M2 = -2

    # Asymptotic behavior of psi near z=0 is: C1*z^Delta_minus + C2*z^Delta_plus
    # Condensation requires the source (C1) to be zero.
    DELTA_PLUS = 2 + np.sqrt(4 + M2)

    # --- 2. Define metric functions for the EGB black hole ---
    def f(z, lambda_gb):
        """Metric function f(z) for EGB black hole (z_h=1, L=1)."""
        sqrt_arg = 1 - 4 * lambda_gb * (1 - z**4)
        if sqrt_arg < 0: return np.nan
        return (1 - np.sqrt(sqrt_arg)) / (2 * lambda_gb)

    def dfdz(z, lambda_gb):
        """Derivative of f(z) with respect to z."""
        # This is a safe way to handle the z=0 and lambda_gb=0 cases
        z = max(z, 1e-12)
        sqrt_arg = 1 - 4 * lambda_gb * (1 - z**4)
        if sqrt_arg <= 0: return np.nan
        return (-4 * z**3) / np.sqrt(sqrt_arg)

    # --- 3. Set up the ODE system for the scalar field psi ---
    def ode_system(z, y, mu):
        psi, dpsi = y
        # Use small epsilons to avoid numerical singularities at the boundaries
        z = max(z, 1e-9)
        z = min(z, 1 - 1e-9)

        f_val = f(z, LAMBDA_GB)
        df_val = dfdz(z, LAMBDA_GB)

        # Coefficients of the ODE: psi'' + P(z)psi' + Q(z)psi = 0
        P_z = (df_val / f_val) - 3 / z
        Q_z = (mu**2 * (1 - z**2)**2) / f_val**2 - M2 / (z**2 * f_val)

        d2psi = -P_z * dpsi - Q_z * psi
        return [dpsi, d2psi]

    # --- 4. Define the shooting function ---
    def shoot(mu):
        z_start = 1 - 1e-5
        z_end = 1e-5

        # Initial conditions from regularity at the horizon: psi'(1) + (M2/4)*psi(1) = 0
        psi_start = 1.0
        dpsi_start = -M2 / 4.0 * psi_start

        sol = solve_ivp(ode_system, [z_start, z_end], [psi_start, dpsi_start], args=(mu,), dense_output=True)

        # Extract solution at the end point near the boundary z=0
        psi_end = sol.sol(z_end)[0]
        dpsi_end = sol.sol(z_end)[1]

        # The zero-source condition corresponds to finding mu where the numerator
        # of the expression for the source coefficient (C1) is zero.
        target_value = z_end * dpsi_end - DELTA_PLUS * psi_end
        return target_value

    # --- 5. Find the critical potential using a root-finding algorithm ---
    try:
        # Search in a reasonable interval based on existing literature
        mu_c = brentq(shoot, 3.0, 5.0)
    except ValueError:
        print("Error: Could not find a root in the specified interval. This may be due to numerical instability or an incorrect search range.")
        return

    # --- 6. Output the result ---
    print(f"For a holographic model with a charged scalar field in a 5D Einstein-Gauss-Bonnet background:")
    print("\n--- Model Parameters ---")
    print(f"Gauss-Bonnet coupling (lambda_GB): {LAMBDA_GB}")
    print(f"Scalar field mass squared (m²): {M2}")
    print("\n--- Result ---")
    print(f"The critical chemical potential for condensation is μ_c = {mu_c:.5f}")

    print("\n--- The Governing Equation ---")
    print("This result is the lowest eigenvalue 'μ' for the differential equation:")
    print("ψ''(z) + (f'(z)/f(z) - 3/z)ψ'(z) + (μ²(1-z²)²/f(z)² - m²/(z²f(z)))ψ(z) = 0")
    print("where the parameters have the values listed above.")


# Run the calculation and print the output
calculate_critical_potential()