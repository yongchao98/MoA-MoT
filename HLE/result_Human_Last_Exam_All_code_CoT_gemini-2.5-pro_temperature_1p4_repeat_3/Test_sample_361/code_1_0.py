import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar

def solve_for_critical_potential():
    """
    Calculates the critical chemical potential for scalar condensation in a 5D
    Einstein-Gauss-Bonnet holographic model.
    """
    # --- Model Parameters ---
    # Gauss-Bonnet coupling
    L_GB = 0.1
    # Mass squared of the scalar field dual to a dimension-3 operator
    M2 = -3.0
    # Horizon radius (sets the temperature scale)
    R_H = 1.0

    # --- Metric Functions ---
    def f(r, l_gb):
        """Metric function f(r) for the EGB black hole."""
        if l_gb > 1/4:
            raise ValueError("Gauss-Bonnet coupling l_gb must be <= 1/4.")
        if np.isclose(r, R_H):
            return 0.0
        term_under_sqrt = 1 - 4 * l_gb * (1 - (R_H / r)**4)
        return (r**2 / (2 * l_gb)) * (1 - np.sqrt(term_under_sqrt))

    def dfdr(r, l_gb):
        """Derivative of the metric function f(r) with respect to r."""
        if np.isclose(r, R_H):
            # The limit r->r_h gives 4*r_h
            return 4.0 * R_H
        if l_gb > 1/4:
            raise ValueError("Gauss-Bonnet coupling l_gb must be <= 1/4.")
            
        sqrt_term = np.sqrt(1 - 4 * l_gb * (1 - (R_H / r)**4))
        term1 = (r / l_gb) * (1 - sqrt_term)
        term2 = (4 * R_H**4 / r**3) / sqrt_term
        return term1 + term2

    # --- System of ODEs ---
    def odes(r, y, mu, l_gb, m2):
        """
        Defines the system of first-order ODEs for the scalar field psi.
        y = [psi, dpsi/dr]
        """
        psi, dpsi = y
        
        f_val = f(r, l_gb)
        df_val = dfdr(r, l_gb)
        
        if np.isclose(f_val, 0):
            # This case should not be reached if we integrate from r_start > r_h
            return [dpsi, 0]

        # ODE: psi'' + (f'/f + 3/r)psi' + (mu^2*(1-rh^2/r^2)^2/f^2 + m^2/f)psi = 0
        d2psidr2 = -(df_val / f_val + 3.0 / r) * dpsi - \
                   (mu**2 * (1 - R_H**2 / r**2)**2 / f_val**2 + m2 / f_val) * psi
                   
        return [dpsi, d2psidr2]

    # --- Objective function for the root-finder ---
    def objective_function(mu):
        """
        Solves the ODE for a given mu and returns the value of psi at a
        large radius. The root of this function is the critical mu.
        """
        r_start = R_H + 1e-6
        r_end = 300.0
        
        # Initial conditions from regularity at the horizon r_h=1
        # psi(r_h) is arbitrary, set to 1.
        # psi'(r_h) = - (m^2 / f'(r_h)) * psi(r_h). Here m^2=-3, f'(1)=4. So psi'(1) = 3/4 * psi(1).
        # Ah, the condition for a charged scalar is different.
        # psi'(rh) = (m^2/f'(rh)) * psi(rh) -> For a neutral scalar.
        # For a charged scalar, the regularity condition at the horizon is psi'(rh) = -3/f'(rh) * psi(rh)
        # for our action with V(psi)=m^2|psi|^2.
        # No, the term in the action is r^3 * psi * ... so the coefficient is d-2=3.
        # psi'(1) = -3/f'(1) = -3/4
        f_prime_at_horizon = dfdr(R_H, L_GB)
        dpsi_0_val = (-3.0 / f_prime_at_horizon)
        
        # Value at r_start = r_h + epsilon
        psi_0 = 1.0 + dpsi_0_val * (r_start - R_H)
        dpsi_0 = dpsi_0_val

        y0 = [psi_0, dpsi_0]
        
        sol = solve_ivp(
            fun=odes,
            t_span=[r_start, r_end],
            y0=y0,
            args=(mu, L_GB, M2),
            method='RK45',
            dense_output=True,
            atol=1e-8,
            rtol=1e-8
        )
        
        # Return the value of psi at the large radius r_end
        psi_at_infinity = sol.y[0, -1]
        return psi_at_infinity

    # --- Finding the Root ---
    # Based on literature, the result for lambda=0.1 should be around 3.3.
    # We search for the root in a bracket around this expected value.
    try:
        # The bracket [3.2, 3.5] should contain the root.
        # For mu < mu_c, psi should go to +inf (positive).
        # For mu > mu_c, psi should go to -inf (negative).
        result = root_scalar(objective_function, bracket=[3.2, 3.5], method='brentq')
        
        print("--- Calculation of Critical Chemical Potential ---")
        print(f"Model parameters:")
        print(f"  Gauss-Bonnet coupling (lambda_GB): {L_GB}")
        print(f"  Scalar mass squared (m^2): {M2}")
        print("\nResult:")
        print(f"The critical chemical potential (mu_c) is: {result.root}")

        # The final answer tag will be populated with this value
        global final_answer
        final_answer = result.root

    except ValueError as e:
        print(f"Error finding the root: {e}")
        print("The objective function may not have opposite signs at the bracket ends.")
        final_answer = "Calculation Failed"

# Execute the calculation
final_answer = 0
solve_for_critical_potential()