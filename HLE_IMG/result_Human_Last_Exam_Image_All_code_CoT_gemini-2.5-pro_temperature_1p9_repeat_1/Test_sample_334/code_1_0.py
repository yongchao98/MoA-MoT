import numpy as np
from scipy.optimize import minimize

def solve_soc_bec():
    """
    Solves the SOC-BEC problem to find the value of n0 * k_R* / k0*.
    """
    
    # Step 1 & 2: Identify parameters based on analysis
    # Base set (Plot n0=3): (delta, Omega, k_R)
    n0 = 3
    delta_0, Omega_0, k_R0 = 2, 4, 1
    
    # Missing parameter set (k_R -> k_R/2 variation)
    delta_star = delta_0
    Omega_star = Omega_0
    k_R_star = k_R0 / 2

    # Step 3: Define the relevant functions for the missing set
    # We need to find k_0 > 0 where d(k*v(k))/dk = 0.
    # This is equivalent to finding the extrema of f(k) = k*v(k).

    def v(k, delta, Omega, k_R):
        """Calculates the group velocity v(k)."""
        # Dispersion relation: E_ = k^2+k_R^2 + d/2 - 0.5*sqrt((4*k*k_R-d)^2+O^2)
        # v(k) = dE_/dk
        denominator = np.sqrt((4 * k * k_R - delta)**2 + Omega**2)
        # The denominator can be zero if k is complex, but for real k it's >= Omega > 0
        if np.isscalar(denominator) and denominator == 0:
            return np.inf 
        
        v_k = 2 * k - (2 * k_R * (4 * k * k_R - delta)) / denominator
        return v_k

    def f_to_minimize(k_array):
        """Function f(k) = k * v(k), to be minimized to find extrema.
           We want to find k_0* > 0 that is an extremum of f(k).
           Numerically, we can find a minimum of -abs(f(k)) to locate non-zero extrema.
           However, finding the root of its derivative is more direct.
        """
        k = k_array[0]
        if k <= 0:  # We are looking for the smallest positive k
            return np.inf
            
        # We need to solve for k where d(k*v(k))/dk = 0.
        # This is k*v'(k) + v(k) = 0.
        
        # Define D = (4*k*k_R - delta)^2 + Omega^2
        D = (4 * k * k_R_star - delta_star)**2 + Omega_star**2
        
        # v(k) as calculated before
        v_k = 2*k - (2*k_R_star * (4*k*k_R_star - delta_star)) / np.sqrt(D)
        
        # v'(k) = d(v(k))/dk
        # v'(k) = 2 - d/dk[ (8*k*k_R^2 - 2*k_R*d) / D^0.5 ]
        # derivative term is (u'w - uw')/w^2
        # u = 8*k*k_R^2 - 2*k_R*d -> u' = 8*k_R^2
        # w = D^0.5 -> w' = 0.5 * D^-0.5 * (2 * (4*k*k_R-d) * 4*k_R) = 4*k_R*(4*k*k_R-d) * D^-0.5
        # ... this simplifies to v'(k) = 2 - (16 * k_R^2 * Omega^2) / D^1.5
        v_prime_k = 2 - (16 * k_R_star**2 * Omega_star**2) / (D**1.5)
        
        # We need to find k that solves k*v'(k) + v(k) = 0
        equation_val = k * v_prime_k + v_k
        
        # The solver minimizes, so we give it the squared value to find a root
        return equation_val**2


    # Initial guess for k > 0
    k_initial_guess = np.array([0.1])
    
    # Use a numerical solver to find the smallest positive k_0*
    # that satisfies d(k*v)/dk = 0.
    result = minimize(f_to_minimize, k_initial_guess, bounds=[(1e-6, None)])
    
    k0_star = result.x[0]

    # Step 4: Final calculation
    final_value = n0 * k_R_star / k0_star

    # Print the equation and the final answer
    print(f"Identified base plot: n_0 = {n0}")
    print(f"Missing parameter set: (delta*, Omega*, k_R*) = ({delta_star}, {Omega_star}, {k_R_star})")
    print(f"Solving for the smallest positive k_0* such that m1(k_0*) + m2(k_0*) = 0.")
    print(f"This is equivalent to solving d(k*v(k))/dk = 0.")
    print(f"Numerically found smallest positive solution: k_0* = {k0_star:.6f}")
    print(f"The final calculation is n0 * k_R* / k0*")
    print(f"{n0} * {k_R_star} / {k0_star:.6f} = {final_value:.6f}")

solve_soc_bec()