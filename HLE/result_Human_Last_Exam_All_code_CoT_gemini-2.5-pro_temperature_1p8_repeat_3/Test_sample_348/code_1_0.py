import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

def main():
    """
    This script calculates the total welfare based on given supply and demand functions.
    """
    
    # Define the demand and supply price functions
    def demand_price(q):
        return 18 * np.exp(-np.arctan(q))

    def supply_price(q):
        # The function is defined for Q^3 > 2, but we are interested in P >= 0.
        # This occurs when Q^3 >= 3.
        # Handle potential domain errors for values outside the integration range.
        if isinstance(q, (np.ndarray, list)):
             # Vectorized operation for quad integrator
            return np.log(np.maximum(q**3 - 2, 1e-9)) # Use a small epsilon to avoid log(0)
        else:
            # Scalar operation for fsolve
            return np.log(q**3 - 2) if q**3 > 2 else -np.inf


    # --- Step 1: Find Market Equilibrium ---
    print("Step 1: Finding Market Equilibrium...")

    # Define the equation to solve for equilibrium: P_demand(Q) - P_supply(Q) = 0
    def equilibrium_equation(q):
        return demand_price(q) - supply_price(q)

    # Solve for the equilibrium quantity Q_E. An initial guess of 4.5 is used based on
    # preliminary analysis showing the intersection is between Q=4 and Q=5.
    q_equilibrium = fsolve(equilibrium_equation, 4.5)[0]

    # Calculate the equilibrium price P_E using the found quantity
    p_equilibrium = demand_price(q_equilibrium)

    print(f"  - Equilibrium Quantity (Q_E): {q_equilibrium:.4f}")
    print(f"  - Equilibrium Price (P_E): {p_equilibrium:.4f}")
    print("-" * 30)

    # --- Step 2: Determine Integration Limits and Calculate Welfare ---
    print("Step 2: Calculating Total Welfare...")
    
    # The market operates for positive prices. For the supply curve, P >= 0 when Q^3 - 2 >= 1, so Q >= 3^(1/3).
    # This sets the lower bound for our integration.
    lower_bound_q = 3**(1/3)
    
    print(f"  - Integration will be performed from Q = {lower_bound_q:.4f} to Q = {q_equilibrium:.4f}")

    # Define the integrand for Consumer Surplus: P_demand(Q) - P_E
    def cs_integrand(q):
        return demand_price(q) - p_equilibrium

    # Define the integrand for Producer Surplus: P_E - P_supply(Q)
    def ps_integrand(q):
        return p_equilibrium - supply_price(q)

    # Calculate Consumer Surplus using numerical integration
    consumer_surplus, _ = quad(cs_integrand, lower_bound_q, q_equilibrium)
    print(f"  - Calculated Consumer Surplus (CS): {consumer_surplus:.4f}")
    
    # Calculate Producer Surplus using numerical integration
    producer_surplus, _ = quad(ps_integrand, lower_bound_q, q_equilibrium)
    print(f"  - Calculated Producer Surplus (PS): {producer_surplus:.4f}")
    print("-" * 30)

    # --- Step 3: Present the Final Result ---
    # Total welfare is the sum of the two surpluses
    total_welfare = consumer_surplus + producer_surplus

    print("Final Result:")
    print("The total welfare is the sum of consumer and producer surplus.")
    print(f"{consumer_surplus:.4f} (CS) + {producer_surplus:.4f} (PS) = {total_welfare:.4f} (Total Welfare)")


if __name__ == "__main__":
    main()
<<<11.1396>>>