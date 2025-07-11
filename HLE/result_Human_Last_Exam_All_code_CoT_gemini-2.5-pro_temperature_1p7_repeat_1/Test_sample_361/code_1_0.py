import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

def solve_holographic_superconductor():
    """
    Calculates the critical chemical potential for a holographic superconductor
    in 5D Einstein-Gauss-Bonnet gravity using the shooting method.
    """

    # Model Parameters
    LAMBDA_GB = 0.1  # Gauss-Bonnet coupling
    M2L2 = -3        # Squared mass of the scalar field (for dual operator dim=3)
    Q = 1            # Charge of the scalar field
    ZH = 1.0         # Horizon position (sets the energy scale)

    # --- Define the background metric functions ---

    def f(z, lambda_gb):
        """EGB blackening factor f(z)."""
        # Ensure the argument of sqrt is non-negative
        sqrt_arg = 1 - 4 * lambda_gb * (1 - (z / ZH)**4)
        if sqrt_arg < 0:
            return -np.inf # Should not happen for lambda_gb <= 1/4
        return (1 / (2 * lambda_gb)) * (1 - np.sqrt(sqrt_arg))

    def df_dz(z, lambda_gb):
        """Derivative of the blackening factor, f'(z)."""
        if z == 0:
            return 0.0
        sqrt_arg = 1 - 4 * lambda_gb * (1 - (z / ZH)**4)
        if sqrt_arg <= 0:
            return -np.inf
        return (4 * z**3) / (ZH**4 * np.sqrt(sqrt_arg))

    # --- Define the ODE system for the scalar field Psi ---

    def ode_system(z, y, mu):
        """
        Defines the system of first-order ODEs.
        y = [Psi, Pi], where Pi = Psi'
        """
        psi, pi = y
        
        # Avoid division by zero at z=0 by handling it as a limit
        if z == 0:
            return [0, 0]

        f_val = f(z, LAMBDA_GB)
        df_val = df_dz(z, LAMBDA_GB)
        
        # Electric potential phi(z)
        phi = mu * (1 - (z/ZH)**2)
        
        # The ODE for Psi(z) converted to a system of 1st order ODEs
        # Psi'' + (f'/f - 3/z)Psi' + (q^2*phi^2/f^2 - m^2/(z^2*f))Psi = 0
        d_psi_dz = pi
        d_pi_dz = -((df_val / f_val) - 3/z) * pi - ((Q**2 * phi**2 / f_val**2) - M2L2 / (z**2 * f_val)) * psi
        
        return [d_psi_dz, d_pi_dz]

    # --- Shooting method objective function ---
    
    def objective_function(mu):
        """
        Solves the ODE for a given mu and returns the value of Psi at the horizon.
        The root of this function is the critical chemical potential mu_c.
        """
        z_start = 1e-6 # A point very close to the boundary z=0
        z_end = ZH - 1e-5 # A point very close to the horizon
        
        # Initial conditions for Psi ~ z^3 near the boundary
        y_start = [z_start**3, 3 * z_start**2]
        
        # Integrate the ODE
        sol = solve_ivp(
            ode_system, 
            [z_start, z_end], 
            y_start, 
            args=(mu,), 
            dense_output=True,
            method='RK45'
        )
        
        # Return the value of Psi at the endpoint near the horizon
        # A solution that diverges will result in a large value here.
        # The root-finder will find where this value crosses zero.
        psi_end = sol.y[0, -1]
        return psi_end

    # --- Find the root to get the critical chemical potential ---

    print("Searching for the critical chemical potential...")
    print("This may take a few seconds.")
    
    # Bracketing the solution. From literature, we expect mu_c to be around 4-5.
    mu_low = 4.5
    mu_high = 5.0
    
    try:
        # Use Brent's method to find the root mu_c where the objective function is zero.
        mu_c = brentq(objective_function, mu_low, mu_high)
        
        print("\n--- Calculation Complete ---")
        print("For the holographic model with the following parameters:")
        print(f"Gauss-Bonnet coupling lambda_GB = {LAMBDA_GB}")
        print(f"Scalar field mass squared m^2*L^2 = {M2L2}")
        print(f"Scalar field charge q = {Q}")
        print("\nThe calculated critical chemical potential is:")
        print(f"mu_c = {mu_c:.8f}")

    except ValueError:
        print(f"Failed to find a root in the interval [{mu_low}, {mu_high}].")
        print("Try adjusting the bracketing interval.")
    except Exception as e:
        print(f"An error occurred during calculation: {e}")

solve_holographic_superconductor()