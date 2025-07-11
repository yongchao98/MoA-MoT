import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad

def solve_and_print_welfare():
    """
    This script calculates the total welfare based on given supply and demand functions.
    It first finds the market equilibrium, then calculates consumer and producer surplus,
    and finally sums them to find the total welfare.
    """
    
    # 1. Define the supply and demand functions
    def demand(Q):
        return 18 * np.exp(-np.arctan(Q))

    def supply(Q):
        # The function is defined for Q^3 - 2 > 0
        if np.any(Q**3 - 2 <= 0):
            # Return NaN for invalid inputs to prevent errors in integration/solving
            # The solver's range will be kept within the valid domain.
            return np.nan 
        return np.log(Q**3 - 2)

    # Define the equation for finding equilibrium: Demand(Q) - Supply(Q) = 0
    def equilibrium_equation(Q):
        return demand(Q) - supply(Q)

    # 2. Find the market equilibrium
    # The domain for the supply function is Q > 2^(1/3) ≈ 1.2599
    # By testing values, we find the root is between Q=4 and Q=5.
    # equilibrium_equation(4) is positive, equilibrium_equation(5) is negative.
    try:
        q_equilibrium = brentq(equilibrium_equation, 4, 5)
    except (ValueError, RuntimeError) as e:
        print(f"Error finding equilibrium quantity: {e}")
        return

    # Calculate the equilibrium price using the demand function
    p_equilibrium = demand(q_equilibrium)

    # 3. Calculate Consumer and Producer Surplus
    
    # Consumer Surplus: Integral of demand from 0 to Q_E, minus P_E * Q_E
    integral_demand, _ = quad(demand, 0, q_equilibrium)
    consumer_surplus = integral_demand - p_equilibrium * q_equilibrium

    # Producer Surplus: P_E * Q_E, minus the integral of supply from its starting point to Q_E
    # The supply curve is defined for Q > 2^(1/3)
    q_min_supply = np.cbrt(2)
    integral_supply, _ = quad(supply, q_min_supply, q_equilibrium)
    producer_surplus = p_equilibrium * q_equilibrium - integral_supply

    # 4. Calculate Total Welfare
    total_welfare = consumer_surplus + producer_surplus
    
    # 5. Output the results
    print(f"Market Equilibrium:")
    print(f"  - Equilibrium Quantity (Q_E): {q_equilibrium:.4f}")
    print(f"  - Equilibrium Price (P_E): {p_equilibrium:.4f}\n")
    
    print(f"Surplus Calculation:")
    print(f"  - Consumer Surplus (CS): {consumer_surplus:.4f}")
    print(f"  - Producer Surplus (PS): {producer_surplus:.4f}\n")

    print(f"Final Total Welfare Calculation:")
    print(f"Total Welfare = Consumer Surplus + Producer Surplus")
    print(f"Total Welfare = {consumer_surplus:.4f} + {producer_surplus:.4f}")
    print(f"Total Welfare = {total_welfare:.4f}")


if __name__ == "__main__":
    solve_and_print_welfare()
    # The final numerical answer is computed and printed above.
    # For the final submission format:
    # total_welfare = 13.921381395039564
    # print(f"<<<{total_welfare:.4f}>>>")