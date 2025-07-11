import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import lambertw
from scipy.optimize import brentq
import warnings

# Suppress warnings that may arise from lambertw function with certain inputs
warnings.filterwarnings("ignore", category=np.ComplexWarning)

def solve_quantum_well_problem():
    """
    This script calculates the energy difference between the first and second
    energy levels of a particle in a specified 3D potential well.
    """
    # 1. Define constants and parameters in convenient units (eV, nm)
    V0_eV = 15.0  # Potential constant in eV
    R_nm = 3.0    # Well radius in nm
    # Pre-calculate hbar^2 / (2m) in eV*nm^2
    HBARC = 197.327      # h-bar * c in eV*nm
    M_E_C2 = 0.510998e6  # Electron rest mass energy in eV
    HBAR2_2M = (HBARC**2) / (2 * M_E_C2)

    # 2. Define the potential energy function V(r)
    def potential_V(r_nm):
        """
        Calculates the potential V(r) in eV at a given radius r in nm.
        """
        # The potential is defined piecewise
        if r_nm < R_nm and r_nm > 0:
            # For 0 < r < R, handle the Lambert W function.
            # The argument to exp(r-R) is assumed to use the numerical value of r and R in nm.
            arg_exp = np.exp(r_nm - R_nm)
            w_val = lambertw(arg_exp).real  # Use the real part of the principal branch
            v_squared = V0_eV + w_val
        elif r_nm >= R_nm:
            # For r >= R
            v_squared = V0_eV * (1.0 - (r_nm / R_nm)**(-2.0))
        else: # r <= 0
            # Handle r=0 case based on the r<R formula
            arg_exp = np.exp(0 - R_nm)
            w_val = lambertw(arg_exp).real
            v_squared = V0_eV + w_val

        # V(r) is the square root of V^2(r). Return a large number if V^2 is negative.
        return np.sqrt(v_squared) if v_squared >= 0 else 1e9

    # 3. Define the ODE system from the radial Schrödinger equation
    def schrodinger_ode(r, y, E_eV):
        """
        Defines the ODE system for the l=0 radial Schrodinger equation.
        y = [u, u'], where u is the radial wavefunction part.
        E_eV is the energy guess in eV.
        """
        u, v = y # y[0] = u(r), y[1] = u'(r)
        du_dr = v
        dv_dr = -(1.0 / HBAR2_2M) * (E_eV - potential_V(r)) * u
        return [du_dr, dv_dr]

    # 4. Define the shooting function
    def wave_function_at_rmax(E_eV):
        """
        Solves the ODE for a given energy E and returns the value of the
        wavefunction u at a large radius, r_max.
        The roots of this function are the energy eigenvalues.
        """
        # Integration range
        r_start = 1e-9  # Start integration slightly away from r=0
        r_max = 8 * R_nm # Choose a sufficiently large radius for the boundary condition

        # Initial conditions: u(r) ~ r near origin, so u'(0) is a constant.
        y0 = [r_start, 1.0]

        # Solve the ODE
        sol = solve_ivp(
            fun=schrodinger_ode,
            t_span=[r_start, r_max],
            y0=y0,
            args=(E_eV,),
            method='RK45'
        )
        
        # Return the final value of u(r)
        return sol.y[0][-1]

    # 5. Find the first two energy eigenvalues (E1, E2) using a root finder
    try:
        # Potential at infinity determines the upper bound for bound state energies
        V_inf = np.sqrt(V0_eV)

        # Find E1 (ground state). Based on analysis, it lies between 0 and V_inf.
        # We provide a reasonable bracket to the root finder.
        E1_bracket = [0.1, 1.5]
        E1 = brentq(wave_function_at_rmax, E1_bracket[0], E1_bracket[1])

        # Find E2 (first excited state). It must be > E1.
        E2_bracket = [E1 + 0.1, V_inf - 0.1]
        E2 = brentq(wave_function_at_rmax, E2_bracket[0], E2_bracket[1])
        
        # 6. Calculate the energy difference
        delta_E = E2 - E1

        # Print the final result including the full equation
        print("Calculation Results:")
        print(f"First energy level (ground state): E1 = {E1:.4f} eV")
        print(f"Second energy level (first excited state): E2 = {E2:.4f} eV")
        print("\nFinal energy difference calculation:")
        print(f"ΔE = {E2:.4f} eV - {E1:.4f} eV = {delta_E:.4f} eV")
        
        global final_answer
        final_answer = delta_E

    except ValueError as e:
        print("Root finding failed. The energy brackets might be incorrect.")
        print(f"Error: {e}")
        final_answer = 0
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        final_answer = 0


# Execute the solver function
final_answer = 0
solve_quantum_well_problem()
print(f"\n<<<{final_answer:.4f}>>>")
