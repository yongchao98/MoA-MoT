import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

def solve_density_profile():
    """
    This function solves for the density profile of a non-ideal gas mixture
    in a gravitational field.

    The method involves the following steps:
    1.  The system is modeled using the principle of local thermodynamic equilibrium. In this framework,
        the chemical potential (mu_i) plus the gravitational potential for each species must be
        constant with height (z). This gives the relation: d(mu_i)/dz = -M_i * g, where M_i is
        the molar mass of species i and g is the gravitational acceleration.
    2.  The chemical potential expression for a van der Waals fluid mixture is used. This makes
        d(mu_i)/dz a function of the local molar concentrations, n_A(z) and n_B(z), and their
        gradients. This relationship forms a system of coupled ordinary differential equations
        (ODEs) for n_A(z) and n_B(z).
    3.  The total number of particles for each gas is known. This constraint serves as a boundary
        condition for the ODE system. The problem is solved as a boundary value problem using
        the shooting method.
    4.  An optimization routine (`scipy.optimize.root`) is employed to find the correct concentrations
        at the bottom of the container (z=0) that satisfy the total particle number constraints.
        It iteratively guesses the initial conditions, solves the ODEs, and adjusts the guess
        to minimize the error.
    5.  With the correct initial conditions determined, the ODEs are solved one final time to obtain
        the precise concentration profiles, n_A(z) and n_B(z).
    6.  Finally, the mass density profile is calculated from the concentration profiles using the
        formula: rho(z) = n_A(z) * M_A + n_B(z) * M_B.
    """
    # --- Parameters ---
    # Gas A
    M_A = 28.0  # g/mol
    N_A = 2.0e23 # number of particles
    a_AA = 2.5   # Pa * m^6 * mol^-2
    b_AA = 0.04  # m^3 * mol^-1

    # Gas B
    M_B = 44.0  # g/mol
    N_B = 1.5e23 # number of particles
    a_BB = 3.6   # Pa * m^6 * mol^-2
    b_BB = 0.05  # m^3 * mol^-1
    
    # Interaction
    a_AB = 3.0   # Pa * m^6 * mol^-2

    # Container and Environment
    A = 0.1      # m^2
    H = 10.0     # m
    T = 500.0    # K
    g = 9.81     # m/s^2

    # --- Constants ---
    N_avo = 6.022e23 # Avogadro's number, mol^-1
    R = 8.314        # Ideal gas constant, J * mol^-1 * K^-1
    
    # --- Derived Parameters in SI units ---
    M_A_kg = M_A / 1000.0 # kg/mol
    M_B_kg = M_B / 1000.0 # kg/mol
    
    N_A_mol = N_A / N_avo # total moles of A
    N_B_mol = N_B / N_avo # total moles of B

    # --- ODE System Definition ---
    # The ODE system is d(y)/dz = J^-1 * g_vec
    # where y = [n_A, n_B], J is the Jacobian of chemical potentials,
    # and g_vec = [-M_A*g, -M_B*g].
    def ode_system(z, y):
        n_A, n_B = y

        # Denominator in vdW expressions
        D = 1.0 - n_A * b_AA - n_B * b_BB
        
        # Check for physical validity. Return a large gradient if invalid
        # to steer the solver away from this region.
        if n_A <= 0 or n_B <= 0 or D <= 1e-9:
            return np.array([1e12, 1e12])
        
        n_tot = n_A + n_B
        
        # Jacobian matrix elements: J_ij = d(mu_i)/d(n_j)
        J_AA = (R * T / n_A) + (2 * R * T * b_AA / D) + (R * T * b_AA**2 * n_tot / D**2) - 2 * a_AA
        J_BB = (R * T / n_B) + (2 * R * T * b_BB / D) + (R * T * b_BB**2 * n_tot / D**2) - 2 * a_BB
        J_AB = (R * T * (b_AA + b_BB) / D) + (R * T * b_AA * b_BB * n_tot / D**2) - 2 * a_AB
        J_BA = J_AB

        J = np.array([[J_AA, J_AB], [J_BA, J_BB]])
        
        # Check if Jacobian is singular
        if abs(np.linalg.det(J)) < 1e-12:
             return np.array([1e12, 1e12])

        J_inv = np.linalg.inv(J)
        g_vec = np.array([-M_A_kg * g, -M_B_kg * g])
        
        return J_inv @ g_vec

    # --- Objective Function for Root Finder (Shooting Method) ---
    def objective_func(y0):
        n_A0, n_B0 = y0
        
        # Check initial guess validity
        if n_A0 <= 0 or n_B0 <= 0 or (1.0 - n_A0 * b_AA - n_B0 * b_BB <= 0):
            return np.array([1e6, 1e6])

        # Solve the initial value problem
        sol = solve_ivp(ode_system, [0, H], y0, method='Radau', t_eval=np.linspace(0, H, 101))
        
        if sol.status != 0: # Solver failed
            return np.array([1e6, 1e6])

        # Integrate profiles to get calculated total moles
        N_A_calc = np.trapz(sol.y[0], sol.t) * A
        N_B_calc = np.trapz(sol.y[1], sol.t) * A

        # Return the difference between calculated and actual total moles
        return np.array([N_A_calc - N_A_mol, N_B_calc - N_B_mol])

    # --- Main Execution ---
    # Initial guess for concentrations at z=0 (n_A0, n_B0)
    # Start with values slightly higher than the average concentration
    avg_n_A = N_A_mol / (A * H)
    avg_n_B = N_B_mol / (A * H)
    initial_guess = np.array([avg_n_A * 1.2, avg_n_B * 1.5])
    
    # Find the correct initial conditions [n_A(0), n_B(0)] using a root finder
    solution_root = root(objective_func, initial_guess, method='hybr', tol=1e-8)
    
    if not solution_root.success:
        print("Error: The root-finding process failed to converge.")
        print("Message:", solution_root.message)
        return

    y0_correct = solution_root.x
    
    # Solve the ODEs one last time with the correct initial conditions
    # Use dense_output=True to get a continuous solution for evaluation at any point
    final_sol = solve_ivp(ode_system, [0, H], y0_correct, method='Radau', dense_output=True)

    # --- Calculate and Print Results ---
    # Get concentrations at z = 0, H/2, H
    n_A_z0, n_B_z0 = final_sol.sol(0)
    n_A_z_mid, n_B_z_mid = final_sol.sol(H/2)
    n_A_zH, n_B_zH = final_sol.sol(H)

    # Calculate mass densities at these heights
    rho_z0 = n_A_z0 * M_A_kg + n_B_z0 * M_B_kg
    rho_z_mid = n_A_z_mid * M_A_kg + n_B_z_mid * M_B_kg
    rho_zH = n_A_zH * M_A_kg + n_B_zH * M_B_kg

    print("The mass density profile rho(z) has been determined.")
    print("Below are the densities at the bottom, middle, and top of the container.\n")
    
    print(f"Density at z = 0.0 m:")
    print(f"  rho(0) = n_A(0) * M_A + n_B(0) * M_B")
    print(f"  rho(0) = {n_A_z0:.4f} mol/m^3 * {M_A_kg:.4f} kg/mol + {n_B_z0:.4f} mol/m^3 * {M_B_kg:.4f} kg/mol = {rho_z0:.4f} kg/m^3\n")

    print(f"Density at z = 5.0 m:")
    print(f"  rho(5.0) = n_A(5.0) * M_A + n_B(5.0) * M_B")
    print(f"  rho(5.0) = {n_A_z_mid:.4f} mol/m^3 * {M_A_kg:.4f} kg/mol + {n_B_z_mid:.4f} mol/m^3 * {M_B_kg:.4f} kg/mol = {rho_z_mid:.4f} kg/m^3\n")

    print(f"Density at z = 10.0 m:")
    print(f"  rho(10.0) = n_A(10.0) * M_A + n_B(10.0) * M_B")
    print(f"  rho(10.0) = {n_A_zH:.4f} mol/m^3 * {M_A_kg:.4f} kg/mol + {n_B_zH:.4f} mol/m^3 * {M_B_kg:.4f} kg/mol = {rho_zH:.4f} kg/m^3")


# Execute the main function when the script is run
if __name__ == '__main__':
    solve_density_profile()