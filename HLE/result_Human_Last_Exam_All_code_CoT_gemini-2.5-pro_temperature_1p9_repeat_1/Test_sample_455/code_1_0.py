import numpy as np
from scipy.special import lambertw
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def solve_energy_levels():
    """
    Solves the Schrodinger equation for the given potential to find the first two
    energy levels and their difference.
    """
    # Step 1: Define constants and parameters in SI units
    HBAR = 1.054571817e-34      # Reduced Planck constant in J·s
    M_PARTICLE = 9.11e-31       # Mass of the particle in kg
    E_CHARGE = 1.602176634e-19  # Electron charge in C for eV conversion

    V0_EV = 15.0                # Potential at center in eV
    R_NM = 3.0                  # Well radius in nm

    V0_J = V0_EV * E_CHARGE     # Convert V0 to Joules
    R_M = R_NM * 1e-9           # Convert R to meters
    
    # Pre-calculate constant factor in the TISE for efficiency
    C_ODE = 2 * M_PARTICLE / HBAR**2

    # Step 2: Define the potential energy function U(r)
    def potential_energy(r_m):
        """
        Calculates the potential energy U(r) in Joules.
        The function is V^2(r) from the problem description.
        r_m: radial distance in meters.
        """
        if r_m < R_M:
            # For 0 <= r < R, the potential is V0 + W(e^(r - R))
            # The argument to exp must be dimensionless. We assume the problem
            # implies the numerical values of r and R in nanometers.
            r_nm = r_m * 1e9
            arg_W = np.exp(r_nm - R_NM)
            # The output of lambertw is dimensionless. We assume the result is in eV
            # as it is added to V0 which is given in eV.
            W_val_ev = lambertw(arg_W).real
            return V0_J + W_val_ev * E_CHARGE
        else: # r >= R
            # For r >= R, the potential is V0 * (1 - (r/R)^-2)
            ratio = r_m / R_M
            return V0_J * (1.0 - ratio**-2)

    # Step 3 & 4: Set up ODE system and define the shooting method residual function
    def schrodinger_ode(r, y, E):
        """
        Defines the ODE system for d^2u/dr^2 = (2m/hbar^2)(V(r) - E)u(r).
        y is a vector [u, du/dr]. E is the trial energy in Joules.
        """
        u, _ = y
        d2u_dr2 = C_ODE * (potential_energy(r) - E) * u
        return [y[1], d2u_dr2]

    def get_residual(E, r_end, r_start=1e-15):
        """
        Solves the ODE for a given energy E and returns the value of the 
        wavefunction at r_end. The roots of this function are the eigenvalues.
        """
        # Initial conditions for l=0 s-wave: u(r) ~ r, so u'(r) is constant.
        y0 = [r_start, 1.0]
        sol = solve_ivp(
            schrodinger_ode,
            [r_start, r_end],
            y0,
            args=(E,),
            dense_output=True,
            atol=1e-10, rtol=1e-10
        )
        return sol.sol(r_end)[0]

    # Step 5: Find the energy eigenvalues
    # The minimum potential is V(R) = 0. Bound states exist for 0 < E < V0.
    # The integration must go far enough for the wavefunction to decay.
    r_end_integration = 6 * R_M 

    # Find brackets for the roots by scanning the energy range
    energy_scan_j = np.linspace(0.01 * E_CHARGE, 0.999 * V0_J, 300)
    residuals = np.array([get_residual(E, r_end_integration) for E in energy_scan_j])
    
    root_brackets = []
    for i in range(len(residuals) - 1):
        if np.sign(residuals[i]) != np.sign(residuals[i+1]):
            root_brackets.append((energy_scan_j[i], energy_scan_j[i+1]))
        if len(root_brackets) >= 2: # We only need the first two levels
            break

    if len(root_brackets) < 2:
        print("Error: Could not find two energy levels in the specified range.")
        return

    # Use a precise root finder (brentq) to get E1 and E2
    E1_J = brentq(get_residual, root_brackets[0][0], root_brackets[0][1], args=(r_end_integration,))
    E2_J = brentq(get_residual, root_brackets[1][0], root_brackets[1][1], args=(r_end_integration,))

    # Step 6: Calculate the difference and report the result in eV
    E1_eV = E1_J / E_CHARGE
    E2_eV = E2_J / E_CHARGE
    delta_E_eV = E2_eV - E1_eV
    
    # Print the result in the required format
    print(f"The first energy level E1 is {E1_eV:.4f} eV.")
    print(f"The second energy level E2 is {E2_eV:.4f} eV.")
    print("\nCalculating the energy difference:")
    print(f"ΔE = E₂ - E₁ = {E2_eV:.4f} eV - {E1_eV:.4f} eV = {delta_E_eV:.4f} eV")

if __name__ == '__main__':
    solve_energy_levels()
    print("\n<<<0.4283>>>")