import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import warnings

def solve_holographic_model():
    """
    This function calculates the critical chemical potential for scalar condensation
    in a 5D Einstein-Gauss-Bonnet holographic model using the shooting method.
    """
    # Suppress integration warnings for cleaner output
    warnings.filterwarnings("ignore", category=UserWarning)

    # --- Step 1: Define Model Parameters ---
    LAMBDA_GB = 0.1  # Gauss-Bonnet coupling
    M2 = -3.0        # Squared mass of the scalar field (m^2 * L^2), for a dual operator of dimension Delta=3
    Q = 1.0          # Charge of the scalar field
    RH = 1.0         # Horizon radius (sets the energy scale, can be set to 1)

    # --- Step 2: Define the Metric Function N(r) ---
    def n_func(r):
        """
        Metric function N(r) for the 5D EGB-AdS black hole.
        We use a Taylor expansion near the horizon to handle numerical precision.
        """
        if np.isclose(r, RH, atol=1e-8):
            # For r -> RH, N(r) ≈ 4 * RH * (r - RH).
            return 4.0 * RH * (r - RH)
        
        arg = 1.0 - 4.0 * LAMBDA_GB * (1.0 - (RH**4 / r**4))
        if arg < 0:
            return np.nan
        return (r**2 / (2.0 * LAMBDA_GB)) * (1.0 - np.sqrt(arg))

    # --- Step 3: Define the ODE System ---
    def odesys(r, y, mu):
        """
        Defines the system of first-order ODEs for the scalar field Ψ.
        y = [Ψ, dΨ/dr]
        """
        psi, dpsi = y
        
        # We use a numerical derivative for N'(r) for robustness.
        eps = 1e-7
        N_r = n_func(r)
        if np.isclose(N_r, 0) or np.isnan(N_r):
             return [dpsi, 0] # Avoid division by zero
        Np_r = (n_func(r + eps) - n_func(r - eps)) / (2*eps)

        # Electric potential Φ(r)
        phi_r = mu * (1.0 - RH**2 / r**2)
        
        # Coefficients in the ODE for Ψ''(r)
        coeff_dpsi = 3.0/r + Np_r/N_r
        coeff_psi = (Q*phi_r/N_r)**2 - M2/N_r
        
        d2psi = -coeff_dpsi * dpsi - coeff_psi * psi
        return [dpsi, d2psi]

    # --- Step 4: Define the Objective Function for Shooting ---
    def objective_func(mu):
        """
        This function returns a value that is zero when the boundary condition at infinity
        is met for a given chemical potential μ.
        """
        if mu <= 0: return 1.0

        r_start = RH + 1e-6
        r_end = 200.0  # A sufficiently large radius to approximate infinity
        r_span = [r_start, r_end]

        # Regularity at the horizon fixes the initial conditions.
        # N'(rH) is analytically known to be 4 for rH=1.
        Np_at_RH = 4.0
        psi0 = 1.0  # Initial value is arbitrary due to linearity
        dpsi0 = M2 / (RH * Np_at_RH)
        
        # Start at r_start using a Taylor expansion
        y0 = [psi0 + dpsi0*(r_start - RH), dpsi0]

        # Solve the initial value problem
        sol = solve_ivp(odesys, r_span, y0, args=(mu,), dense_output=True, method='RK45')
        
        if sol.status != 0: return 1e12 # Penalize failed integration

        # At large r, Ψ(r) ≈ C₁/r + C₃/r³. Condensation requires C₁=0.
        # The target function is constructed to be proportional to C₁.
        y_end = sol.sol(r_end)
        psi_end, dpsi_end = y_end[0], y_end[1]
        target_val = r_end**3 * dpsi_end + 3.0 * r_end**2 * psi_end
        return target_val

    # --- Step 5: Solve and Print the Result ---
    print("In the context of bottom-up holographic models based on a D3/D7 configuration,")
    print("the critical chemical potential for scalar condensation is determined by solving")
    print("the equation of motion for the dual scalar field Ψ in the EGB-AdS background:")
    print("\n  Ψ''(r) + (3/r + N'/N)Ψ' + (q²Φ²/N² - m²/N)Ψ = 0\n")
    print("The solution is found numerically for the following parameters:")
    
    equation_params = {
        "Gauss-Bonnet coupling (λ_GB)": LAMBDA_GB,
        "Scalar field mass squared (m²L²)": M2,
        "Scalar field charge (q)": Q
    }
    for param, value in equation_params.items():
        print(f"  - {param:<30}: {value}")
        
    try:
        # Find the root μ_c within a plausible bracket [1.0, 1.5]
        mu_c = brentq(objective_func, 1.0, 1.5, xtol=1e-6, rtol=1e-6)
        print("\nThe numerical solution for the critical chemical potential is:")
        print(f"  μ_c = {mu_c:.4f}")
        return mu_c
    except ValueError:
        print("\nCould not find the root within the specified bracket.")
        return None

if __name__ == '__main__':
    solve_holographic_model()