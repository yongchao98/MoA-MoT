import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

def main():
    """
    This script calculates the total welfare given supply and demand functions.
    """
    # Step 1: Define the supply and demand functions
    def demand_P(Q):
        """Demand function: P = 18 * e^(-arctan(Q))"""
        return 18 * np.exp(-np.arctan(Q))

    def supply_P(Q):
        """Supply function: P = ln(Q^3 - 2)"""
        # The function is defined for Q^3 - 2 > 0
        with np.errstate(all='ignore'):
            return np.log(Q**3 - 2)

    def equilibrium_equation(Q):
        """Equation to solve for equilibrium: Demand P - Supply P = 0"""
        # This check helps the solver avoid the undefined region of the supply function.
        if Q**3 - 2 <= 0:
            return np.inf
        return demand_P(Q) - supply_P(Q)

    # Step 2: Find the market equilibrium numerically
    # An initial guess of Q=3 is reasonable from plotting or testing values.
    initial_guess_Q = 3.0
    # Use fsolve to find the equilibrium quantity (Q_e)
    equilibrium_quantity, = fsolve(equilibrium_equation, initial_guess_Q)
    
    # Calculate the equilibrium price (P_e)
    equilibrium_price = demand_P(equilibrium_quantity)

    # Step 3: Calculate the total welfare via numerical integration
    # The integral of the demand function from 0 to Q_e
    integral_demand, _ = quad(demand_P, 0, equilibrium_quantity)

    # The integral of the supply function from its starting point (Q^3 - 2 > 0) to Q_e
    supply_lower_bound = 2**(1/3)
    integral_supply, _ = quad(supply_P, supply_lower_bound, equilibrium_quantity)
    
    # Total welfare is the area between the curves
    total_welfare = integral_demand - integral_supply
    
    # Step 4: Print the results in a clear format
    print("--- Market Equilibrium and Total Welfare Calculation ---")
    print(f"\n1. Solved for market equilibrium where Supply equals Demand:")
    print(f"   ln(Q_e^3 - 2) = 18 * e^(-arctan(Q_e))")
    print(f"\n   Equilibrium Quantity (Q_e) = {equilibrium_quantity:.4f}")
    print(f"   Equilibrium Price (P_e) = {equilibrium_price:.4f}")
    
    print("\n2. Calculated Total Welfare by integrating the area between the curves:")
    print(f"   Total Welfare = [Integral of Demand from 0 to Q_e] - [Integral of Supply from 2^(1/3) to Q_e]")
    print("\n   The final welfare equation with calculated values is:")
    print(f"   Total Welfare = {integral_demand:.4f} - {integral_supply:.4f}")
    print(f"   Total Welfare = {total_welfare:.4f}")

if __name__ == "__main__":
    main()
<<<25.8016>>>