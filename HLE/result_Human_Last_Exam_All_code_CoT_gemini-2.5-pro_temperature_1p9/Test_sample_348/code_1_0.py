import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

def main():
    """
    This script calculates the total welfare based on given supply and demand functions.
    """

    # 1. Define supply and demand functions
    def demand_p(q):
        """Calculates the price for a given quantity demanded."""
        return 18 * np.exp(-np.arctan(q))

    def supply_p(q):
        """Calculates the price for a given quantity supplied."""
        # The function is defined for Q^3 - 2 > 0.
        # This guard prevents math errors for values outside the domain.
        if np.any(q**3 <= 2):
            return np.nan
        return np.log(q**3 - 2)

    # 2. Find Market Equilibrium
    def equilibrium_difference(q):
        """Difference between demand and supply price. We find the root of this function."""
        return demand_p(q) - supply_p(q)

    # The supply function is defined for Q > 2^(1/3) ~= 1.26.
    # We use an initial guess greater than this, like 2.0, to find the root.
    try:
        q_equilibrium = fsolve(equilibrium_difference, 2.0)[0]
    except Exception as e:
        print(f"Could not find the equilibrium point. Error: {e}")
        return

    p_equilibrium = demand_p(q_equilibrium)

    print(f"Step 1: Market Equilibrium")
    print(f"Equilibrium Quantity (Q_E): {q_equilibrium:.4f}")
    print(f"Equilibrium Price (P_E): {p_equilibrium:.4f}\n")

    # 3. Calculate Consumer and Producer Surplus
    # Consumer Surplus: Integral of demand from 0 to Q_E, minus P_E * Q_E
    integral_demand, _ = quad(demand_p, 0, q_equilibrium)
    consumer_surplus = integral_demand - (p_equilibrium * q_equilibrium)

    # Producer Surplus: P_E * Q_E minus the integral of supply.
    # The supply curve starts where Q^3 - 2 > 0, so Q_min = 2^(1/3)
    q_min_supply = 2**(1/3)
    integral_supply, _ = quad(supply_p, q_min_supply, q_equilibrium)
    producer_surplus = (p_equilibrium * q_equilibrium) - integral_supply

    print(f"Step 2: Surplus Calculation")
    print(f"Consumer Surplus (CS): {consumer_surplus:.4f}")
    print(f"Producer Surplus (PS): {producer_surplus:.4f}\n")

    # 4. Calculate Total Welfare
    total_welfare = consumer_surplus + producer_surplus

    print(f"Step 3: Total Welfare Calculation")
    # Outputting each number in the final equation as requested.
    print(f"Total Welfare = CS + PS")
    print(f"Total Welfare = {consumer_surplus:.4f} + {producer_surplus:.4f}")
    print(f"Total Welfare = {total_welfare:.4f}")


if __name__ == "__main__":
    main()
<<<11.4552>>>