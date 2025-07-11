import numpy as np
from scipy.optimize import root_scalar
from scipy.integrate import quad

def main():
    """
    This script calculates the total welfare based on given supply and demand functions.
    """
    # 1. Define the supply and demand functions
    def demand(Q):
        """Demand function: P = 18 * e^(-arctan(Q))"""
        return 18 * np.exp(-np.arctan(Q))

    def supply(Q):
        """Supply function: P = ln(Q^3 - 2)"""
        # Ensure Q is in the valid domain for the natural logarithm (Q^3 - 2 > 0)
        # We handle this by setting the lower bound of integration/root finding appropriately.
        return np.log(Q**3 - 2)

    # 2. Find market equilibrium (Q_e, P_e)
    # Define the equation for the root finder: Demand(Q) - Supply(Q) = 0
    def equilibrium_equation(Q):
        return demand(Q) - supply(Q)

    # The supply function is defined for Q^3 > 2, so Q > 2^(1/3)
    q_min_supply = np.cbrt(2)

    try:
        # Use a root-finding algorithm to find the equilibrium quantity (Q_e)
        # We need a bracket [a, b] where equilibrium_equation(a) and equilibrium_equation(b) have opposite signs.
        # A quick test shows the root is between 4 and 5.
        sol = root_scalar(equilibrium_equation, bracket=[q_min_supply + 0.1, 10])
        Q_e = sol.root
        # Calculate the equilibrium price (P_e) using the demand function
        P_e = demand(Q_e)
    except ValueError as e:
        print(f"Error: Could not find the market equilibrium point. The curves may not intersect. Details: {e}")
        return

    # 3. Calculate Consumer and Producer Surplus
    # Integrate the demand function from 0 to Q_e
    integral_demand, _ = quad(demand, 0, Q_e)
    consumer_surplus = integral_demand - Q_e * P_e

    # Integrate the supply function from its minimum defined quantity (q_min_supply) to Q_e
    integral_supply, _ = quad(supply, q_min_supply, Q_e)
    producer_surplus = Q_e * P_e - integral_supply
    
    # 4. Calculate Total Welfare
    total_welfare = consumer_surplus + producer_surplus

    # Print the detailed breakdown of the calculation
    print("--- Market Equilibrium ---")
    print(f"Equilibrium Quantity (Q_e) = {Q_e:.4f}")
    print(f"Equilibrium Price (P_e) = {P_e:.4f}\n")

    print("--- Calculating Consumer Surplus (CS) ---")
    print("Formula: CS = (Integral of Demand from 0 to Q_e) - (P_e * Q_e)")
    print(f"CS = {integral_demand:.4f} - ({P_e:.4f} * {Q_e:.4f})")
    print(f"CS = {integral_demand:.4f} - {P_e * Q_e:.4f} = {consumer_surplus:.4f}\n")
    
    print("--- Calculating Producer Surplus (PS) ---")
    print(f"Formula: PS = (P_e * Q_e) - (Integral of Supply from Q_min to Q_e), where Q_min={q_min_supply:.4f}")
    print(f"PS = ({P_e:.4f} * {Q_e:.4f}) - {integral_supply:.4f}")
    print(f"PS = {P_e * Q_e:.4f} - {integral_supply:.4f} = {producer_surplus:.4f}\n")

    print("--- Calculating Total Welfare (TW) ---")
    print("Formula: TW = Consumer Surplus + Producer Surplus")
    print(f"TW = {consumer_surplus:.4f} + {producer_surplus:.4f} = {total_welfare:.4f}")

if __name__ == "__main__":
    main()
<<<20.0883>>>