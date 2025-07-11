import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

def main():
    """
    This script calculates the total welfare at market equilibrium given
    supply and demand functions.
    """

    # Define the demand function P = 18 * exp(-arctan(Q))
    def demand_p(q):
        return 18 * np.exp(-np.arctan(q))

    # Define the supply function P = ln(Q^3 - 2)
    def supply_p(q):
        # The function is defined for Q > cbrt(2). Return np.inf for invalid inputs
        # to guide the solver.
        if q**3 - 2 <= 0:
            return np.inf
        return np.log(q**3 - 2)

    # Define the equation to find equilibrium: Supply_P(Q) - Demand_P(Q) = 0
    def equilibrium_equation(q):
        return supply_p(q) - demand_p(q)

    # --- Step 1: Find Market Equilibrium ---
    # We need an initial guess for the solver. Let's try a value like 4.
    initial_guess = 4.0
    # Use fsolve to find the equilibrium quantity Q_E
    try:
        Q_E = fsolve(equilibrium_equation, initial_guess)[0]
    except Exception as e:
        print(f"Could not find equilibrium point. Error: {e}")
        return

    # Calculate the equilibrium price P_E using the demand function
    P_E = demand_p(Q_E)
    
    # --- Step 2: Calculate Consumer Surplus (CS) ---
    # CS = Integral(Demand) from 0 to Q_E - (P_E * Q_E)
    integral_demand, _ = quad(demand_p, 0, Q_E)
    consumer_surplus = integral_demand - (P_E * Q_E)

    # --- Step 3: Calculate Producer Surplus (PS) ---
    # PS = (P_E * Q_E) - Integral(Supply) from Q_min to Q_E
    # The economically relevant supply starts at P=0, which is when Q^3 - 2 = 1, so Q = cbrt(3).
    q_min_supply = np.cbrt(3)
    
    integral_supply, _ = quad(supply_p, q_min_supply, Q_E)
    producer_surplus = (P_E * Q_E) - integral_supply

    # --- Step 4: Calculate Total Welfare ---
    total_welfare = consumer_surplus + producer_surplus

    # --- Print Results ---
    print("Calculation of Total Welfare")
    print("=" * 40)
    print("1. Market Equilibrium:")
    print(f"   Equilibrium Quantity (Q_E) = {Q_E:.4f}")
    print(f"   Equilibrium Price (P_E)    = {P_E:.4f}")
    print("-" * 40)

    print("2. Consumer Surplus (CS):")
    print(f"   CS = Integral from 0 to {Q_E:.4f} of P_D(Q) dQ - ({P_E:.4f} * {Q_E:.4f})")
    print(f"   CS = {integral_demand:.4f} - {P_E * Q_E:.4f}")
    print(f"   CS = {consumer_surplus:.4f}")
    print("-" * 40)
    
    print("3. Producer Surplus (PS):")
    print(f"   PS = ({P_E:.4f} * {Q_E:.4f}) - Integral from {q_min_supply:.4f} to {Q_E:.4f} of P_S(Q) dQ")
    print(f"   PS = {P_E * Q_E:.4f} - {integral_supply:.4f}")
    print(f"   PS = {producer_surplus:.4f}")
    print("-" * 40)

    print("4. Total Welfare (TW):")
    print(f"   TW = CS + PS")
    print(f"   TW = {consumer_surplus:.4f} + {producer_surplus:.4f}")
    print(f"   TW = {total_welfare:.4f}")
    print("=" * 40)

    # The final answer in the requested format
    # The problem did not specify the required precision, so I will use 4 decimal places.
    print(f"\nFinal Answer (Total Welfare): {total_welfare:.4f}")


if __name__ == "__main__":
    main()
