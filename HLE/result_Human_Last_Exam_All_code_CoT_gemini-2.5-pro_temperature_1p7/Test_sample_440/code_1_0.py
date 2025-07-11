import numpy as np
from scipy.integrate import solve_ivp, trapezoid
from scipy.optimize import root
import sys

def solve_density_profile():
    """
    Calculates the density profile of a two-component non-ideal gas mixture in a gravitational field.
    The primary result returned is the mass density at the bottom of the container, rho(0).
    """
    # Part 1: Constants and Parameters
    # ------------------------------------
    # Physical constants
    N_A_const = 6.02214076e23  # Avogadro's number, mol^-1
    R = 8.314462618           # Ideal gas constant, J/(mol·K)
    k_B = R / N_A_const       # Boltzmann constant, J/K

    # System parameters from the problem
    A = 0.1                   # Cross-sectional area, m^2
    H = 10.0                  # Height, m
    T = 500.0                 # Temperature, K
    g = 9.81                  # Gravitational acceleration, m/s^2
    N_A_target = 2.0e23       # Total particles of Gas A
    N_B_target = 1.5e23       # Total particles of Gas B

    # --- Gas Properties (converted to per-particle SI units) ---
    # Gas A
    M_A_molar = 28e-3         # kg/mol
    m_A = M_A_molar / N_A_const  # kg/particle
    a_AA_molar = 2.5          # Pa·m^6·mol^-2
    b_AA_molar = 0.04         # m^3·mol^-1
    a_AA = a_AA_molar / N_A_const**2
    b_AA = b_AA_molar / N_A_const

    # Gas B
    M_B_molar = 44e-3         # kg/mol
    m_B = M_B_molar / N_A_const  # kg/particle
    a_BB_molar = 3.6          # Pa·m^6·mol^-2
    b_BB_molar = 0.05         # m^3·mol^-1
    a_BB = a_BB_molar / N_A_const**2
    b_BB = b_BB_molar / N_A_const

    # Interaction A-B
    a_AB_molar = 3.0          # Pa·m^6·mol^-2
    a_AB = a_AB_molar / N_A_const**2
    
    # Pack parameters for easy access in functions
    params = locals()

    # Part 2: ODE System Definition
    # -------------------------------
    def ode_system(z, y, p):
        nA, nB = y
        
        # Check for unphysical densities (volume fraction >= 1 or negative densities)
        b_term_sum = nA * p['b_AA'] + nB * p['b_BB']
        if b_term_sum >= 1.0 or nA <= 0 or nB <= 0:
            return [np.inf, np.inf] # Return a large gradient to steer the solver away
            
        f = 1.0 - b_term_sum
        f2 = f**2
        
        # Matrix M = d(mu_i)/d(n_k), where mu is the chemical potential
        M11 = (p['k_B']*p['T']/nA) + (p['k_B']*p['T']*p['b_AA']/f) + (p['k_B']*p['T']*p['b_AA']**2/f2) - 2*p['a_AA']
        M12 = (p['k_B']*p['T']*p['b_BB']/f) + (p['k_B']*p['T']*p['b_AA']*p['b_BB']/f2) - 2*p['a_AB']
        M21 = (p['k_B']*p['T']*p['b_AA']/f) + (p['k_B']*p['T']*p['b_BB']*p['b_AA']/f2) - 2*p['a_AB']
        M22 = (p['k_B']*p['T']/nB) + (p['k_B']*p['T']*p['b_BB']/f) + (p['k_B']*p['T']*p['b_BB']**2/f2) - 2*p['a_BB']
        
        detM = M11 * M22 - M12 * M21
        if abs(detM) < 1e-150: # Avoid numerical singularity
            return [np.inf, np.inf]

        # Inverse of M
        Minv11, Minv12, Minv21, Minv22 = M22/detM, -M12/detM, -M21/detM, M11/detM
        
        # Gravity vector on the RHS of the ODE system
        g_vec = np.array([-p['m_A'] * p['g'], -p['m_B'] * p['g']])
        
        # Derivatives dn/dz = M^-1 * (-m*g)
        dnA_dz = Minv11 * g_vec[0] + Minv12 * g_vec[1]
        dnB_dz = Minv21 * g_vec[0] + Minv22 * g_vec[1]
        
        return [dnA_dz, dnB_dz]

    # Part 3: Objective Function for Root Finder
    # ------------------------------------------
    def objective_function(n0, p):
        nA_0, nB_0 = n0
        # Physicality check for initial guess from the root finder
        if nA_0 <= 0 or nB_0 <= 0 or (nA_0 * p['b_AA'] + nB_0 * p['b_BB'] >= 1.0):
            return [1e10, 1e10] # Return large error

        sol = solve_ivp(
            fun=ode_system, t_span=[0, p['H']], y0=n0, args=(p,),
            method='RK45', dense_output=True, rtol=1e-7, atol=1e-10)

        if sol.status != 0: return [1e10, 1e10] # Solver failed

        z_points = np.linspace(0, p['H'], 200) # Integration points
        n_profiles = sol.sol(z_points)
        nA_profile, nB_profile = n_profiles[0], n_profiles[1]

        if np.any(nA_profile <= 0) or np.any(nB_profile <= 0): return [1e10, 1e10]

        N_A_calc = trapezoid(nA_profile, z_points) * p['A']
        N_B_calc = trapezoid(nB_profile, z_points) * p['A']
        
        # Return the difference between calculated and target particle numbers
        return [N_A_calc - p['N_A_target'], N_B_calc - p['N_B_target']]

    # Part 4: Solve for Initial Conditions
    # ------------------------------------
    V_total = A * H
    nA_avg = N_A_target / V_total
    nB_avg = N_B_target / V_total
    initial_guess = [nA_avg, nB_avg] 

    solution = root(objective_function, initial_guess, args=(params,), method='hybr', tol=1e-9)

    if not solution.success:
        print("Error: The root-finding algorithm failed to converge.", file=sys.stderr)
        print(f"Message from solver: {solution.message}", file=sys.stderr)
        sys.exit(1)
        
    nA_0, nB_0 = solution.x

    # Part 5: Final Calculation and Output
    # ------------------------------------
    rho_0 = nA_0 * m_A + nB_0 * m_B

    print("The mass density at the bottom of the container, rho(z=0), is calculated from the equilibrium number densities at z=0.")
    print("The governing equation is: rho(0) = n_A(0) * m_A + n_B(0) * m_B")
    print("\nCalculated values:")
    print(f"Equilibrium number density of Gas A at z=0: n_A(0) = {nA_0:.4e} particles/m^3")
    print(f"Mass of a Gas A particle: m_A = {m_A:.4e} kg")
    print(f"Equilibrium number density of Gas B at z=0: n_B(0) = {nB_0:.4e} particles/m^3")
    print(f"Mass of a Gas B particle: m_B = {m_B:.4e} kg")
    print("\nFinal Result:")
    print(f"rho(0) = ({nA_0:.4e}) * ({m_A:.4e}) + ({nB_0:.4e}) * ({m_B:.4e})")
    print(f"rho(0) = {rho_0:.6f} kg/m^3")
    print(f"<<<{rho_0:.6f}>>>")

if __name__ == "__main__":
    solve_density_profile()