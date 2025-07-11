import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar

def calculate_critical_potential():
    """
    Calculates the critical chemical potential for scalar condensation in a D3/D7
    holographic model with an Einstein-Gauss-Bonnet background.
    """
    # 1. Define model parameters
    LAMBDA_GB = 0.1
    # We are looking for the proportionality constant k in mu_c = k * sqrt(rho).
    # We can find k by setting rho=1 and finding the corresponding mu.
    RHO = 1.0
    M2 = -3.0  # Mass-squared for the D7-brane embedding scalar (m^2*L^2 = -3)

    # 2. Calculate derived constants
    c_GB_sq = (1 - np.sqrt(1 - 4 * LAMBDA_GB)) / (2 * LAMBDA_GB)
    c_GB = np.sqrt(c_GB_sq)

    # Exponents of the two independent solutions near the boundary z=0
    # s^2 - 4s + m^2 = 0 => s = 2 +/- sqrt(4-m^2)
    s_source = 2 - np.sqrt(4 - M2)      # Corresponds to the source term
    s_vev = 2 + np.sqrt(4 - M2)         # Corresponds to the VEV term

    # 3. Define the ODE system for the shooting method
    def ode_system(z, Y, mu):
        """
        Defines the system of first-order ODEs. Y = [psi, dpsi/dz].
        Original EOM: psi'' - (3/z)psi' + V(z)psi = 0
        """
        psi, dpsi = Y
        # The potential term V(z) in the EOM
        potential_term = (1 / c_GB_sq) * (M2 - (mu - RHO * z**2)**2)
        # Rearrange EOM to solve for psi'': psi'' = (3/z)psi' - V(z)psi
        d2psi = (3 / z) * dpsi - potential_term * psi
        return [dpsi, d2psi]

    # 4. Define the objective function for the root-finder
    def objective_function(mu):
        """
        Solves the ODE for a given mu and returns a value that should be zero
        when the boundary condition at z=0 is met.
        """
        z_max = 8.0   # A large value of z to start integration
        z_min = 1e-5  # A small value of z to check the boundary condition

        # Asymptotic behavior at large z gives the initial conditions
        # psi ~ exp(-rho * z^2 / (2 * c_GB))
        psi_init = np.exp(-RHO * z_max**2 / (2 * c_GB))
        dpsi_init = (-RHO * z_max / c_GB) * psi_init

        # Numerically integrate the ODE from z_max down to z_min
        sol = solve_ivp(
            ode_system,
            [z_max, z_min],
            [psi_init, dpsi_init],
            args=(mu,),
            dense_output=True,
            method='RK45'
        )

        # Extract the solution at z_min
        psi_final = sol.sol(z_min)[0]
        dpsi_final = sol.sol(z_min)[1]

        # The general solution near z=0 is psi = c1*z^s_source + c2*z^s_vev.
        # We need to find mu such that the coefficient c1 of the source term is zero.
        # c1 is proportional to (s_vev * psi - z * dpsi). We want this to be zero.
        target_value = s_vev * psi_final - z_min * dpsi_final
        return target_value

    # 5. Find the root of the objective function
    # Based on literature (e.g., Pan et al., PRD 81, 106007 (2010)), the value for
    # lambda_GB=0.1 should be around 4.1. We use this to set a search bracket.
    try:
        sol = root_scalar(objective_function, bracket=[3.8, 4.4], method='brentq')
        mu_c_over_sqrt_rho = sol.root
        
        print("Calculation successful.")
        print(f"The model is a D3/D7 system in 5D Einstein-Gauss-Bonnet gravity with λ_GB = {LAMBDA_GB}.")
        print("The critical chemical potential (μ_c) for condensation is related to the charge density (ρ) by:")
        print(f"μ_c = k * sqrt(ρ)")
        print("\nThe calculated value for the proportionality constant 'k' is:")
        print(f"k = {mu_c_over_sqrt_rho:.4f}")
        return mu_c_over_sqrt_rho

    except ValueError:
        print("Error: Could not find the root within the specified bracket.")
        print("The solution might be outside the range [3.8, 4.4] or the model is unstable.")
        return None

if __name__ == '__main__':
    result = calculate_critical_potential()
