import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar

def solve_critical_potential():
    """
    Calculates the critical chemical potential for scalar condensation in a 5D
    Einstein-Gauss-Bonnet holographic model using a numerical shooting method.
    """
    # --- Model Parameters ---
    # Gauss-Bonnet coupling
    LAMBDA_GB = 0.1
    # Scalar field mass squared (for dual operator with dimension Delta=3 in d=4)
    MASS_SQUARED = -3.0
    # AdS radius (set to 1)
    L = 1.0
    # Horizon position (sets the temperature scale, set to 1)
    ZH = 1.0

    # --- EGB Metric Functions ---
    def f(z, lgb, zh):
        """EGB blackening function f(z)."""
        term = 1 - 4 * lgb * (1 - (z / zh)**4)
        if term < 0:
            raise ValueError(f"sqrt argument is negative at z={z}. Check parameters.")
        # The form below is numerically stable for small lgb.
        return (1 - np.sqrt(term)) / (2 * lgb)

    def f_prime(z, lgb, zh):
        """Derivative of f(z) with respect to z."""
        if np.isclose(z, 0):
            return 0.0
        term = 1 - 4 * lgb * (1 - (z / zh)**4)
        if term <= 0:
           raise ValueError(f"sqrt argument is zero or negative in f_prime at z={z}.")
        return (4 * (z / zh)**3 / zh) / np.sqrt(term)

    # --- ODE System for the Scalar Field Psi ---
    def ode_system(z, y, mu):
        """
        Defines the second-order ODE for psi(z) as a system of two first-order ODEs.
        y[0] = psi, y[1] = dpsi/dz.
        """
        psi, dpsi_dz = y
        
        # Metric functions at current z
        f_val = f(z, LAMBDA_GB, ZH)
        fp_val = f_prime(z, LAMBDA_GB, ZH)
        
        # Potential A_t(z) at the critical point
        At = mu * (1 - (z/ZH)**2)

        # ODE coefficients: psi'' + A*psi' + B*psi = 0
        # A = (f'/f - 3/z)
        # B = (A_t^2/f^2 - m^2/(z^2*f))
        A = (fp_val / f_val) - (3 / z)
        B = (At**2 / f_val**2) - (MASS_SQUARED / (z**2 * f_val))
        
        d2psi_dz2 = -A * dpsi_dz - B * psi
        return [dpsi_dz, d2psi_dz2]

    # --- Objective Function for the Root Finder ---
    def objective(mu):
        """
        Solves the ODE for a given mu and returns the value of the non-normalizable
        mode at the boundary. The root of this function is the critical potential mu_c.
        """
        # Integrate from near the horizon to near the boundary
        z_start = ZH - 1e-6
        z_end = 1e-6
        
        # Initial conditions at the horizon (z=z_h) from regularity analysis
        # psi(z_h) = 1 (arbitrary normalization)
        # psi'(z_h) = m^2 * L^2 / 4
        y0 = [1.0, MASS_SQUARED / 4.0]
        
        # Solve the initial value problem
        sol = solve_ivp(
            ode_system, 
            [z_start, z_end], 
            y0, 
            args=(mu,),
            method='RK45'
        )
        
        # The value to be zeroed is the coefficient of the non-normalizable mode z^(d-Delta)=z^1.
        # At small z, psi(z) ~ c1*z + c2*z^3. We want to find mu such that c1=0.
        # psi(z_end) is dominated by the c1*z_end term. So we seek mu for which psi(z_end)=0.
        psi_end = sol.y[0, -1]
        return psi_end

    # --- Find the Critical Potential ---
    # Bracket the root based on known physics (mu_c increases with lambda_GB)
    # The value for lambda_GB=0 is ~2.5, so for 0.1 it should be higher.
    try:
        result = root_scalar(objective, bracket=[3.0, 4.0], method='brentq')
        mu_c = result.root
        print(f"The critical chemical potential (for z_h=1) is: {mu_c}")
        return mu_c
    except ValueError as e:
        print(f"Could not find the root in the specified bracket. Error: {e}")
        print("Please try adjusting the bracket in the root_scalar call.")
        return None

if __name__ == '__main__':
    critical_potential = solve_critical_potential()
    if critical_potential is not None:
        print(f"\nFinal Answer for the critical chemical potential, \u03BCc:")
        print(f"<<<{critical_potential:.10f}>>>")
