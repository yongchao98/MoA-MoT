import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

def main():
    """
    This script calculates the total welfare based on given supply and demand functions.
    """

    # 1. Define Supply and Demand Functions
    def supply_func(Q_S):
        """Supply function P = ln(Q_S^3 - 2)"""
        # The function is defined for Q_S^3 - 2 > 0, i.e., Q_S > 2^(1/3)
        # We need to handle cases where the input is outside this domain.
        if isinstance(Q_S, (np.ndarray, list)):
            # Create a safe copy to avoid modifying the original array
            result = np.full_like(Q_S, np.nan, dtype=np.float64)
            valid_indices = Q_S**3 > 2
            result[valid_indices] = np.log(Q_S[valid_indices]**3 - 2)
            return result
        else: # Scalar input
            if Q_S**3 <= 2:
                return np.nan # Undefined
            return np.log(Q_S**3 - 2)

    def demand_func(Q_D):
        """Demand function P = 18 * e^(-arctan(Q_D))"""
        return 18 * np.exp(-np.arctan(Q_D))

    # 2. Find Market Equilibrium
    def equilibrium_equation(Q_arr):
        """Equation to solve for equilibrium: Demand(Q) - Supply(Q) = 0"""
        Q = Q_arr[0]
        # The supply function is only defined for Q > 2^(1/3).
        # If Q is outside this domain, there's no equilibrium. We return a large
        # value to guide the solver towards the valid domain.
        if Q**3 <= 2:
            return 1e9  # A large residual
        return demand_func(Q) - supply_func(Q)

    # An initial guess for Q* (must be in the valid domain for supply)
    initial_guess_Q = 4.0
    Q_star = fsolve(equilibrium_equation, initial_guess_Q)[0]

    # Calculate P* using the demand function (it's simpler)
    P_star = demand_func(Q_star)

    # 3. Calculate Welfare
    
    # Calculate Consumer Surplus (CS)
    # CS = Integral(Demand) from 0 to Q* - P* * Q*
    integral_demand, _ = quad(demand_func, 0, Q_star)
    consumer_surplus = integral_demand - (P_star * Q_star)

    # Calculate Producer Surplus (PS)
    # PS = P* * Q* - Integral(Supply) from its Q-intercept to Q*
    # The supply curve is defined for Q > 2^(1/3). This is its effective start.
    supply_q_intercept = (2)**(1/3)
    
    # We only integrate where the supply curve exists
    integral_supply, _ = quad(supply_func, supply_q_intercept, Q_star)
    producer_surplus = (P_star * Q_star) - integral_supply

    # Calculate Total Welfare
    total_welfare = consumer_surplus + producer_surplus

    # 4. Output the results
    print("--- Market Equilibrium ---")
    print(f"Equilibrium Quantity (Q*): {Q_star:.4f}")
    print(f"Equilibrium Price (P*):    {P_star:.4f}\n")
    
    print("--- Welfare Calculation ---")
    print(f"Consumer Surplus: {consumer_surplus:.4f}")
    print(f"Producer Surplus: {producer_surplus:.4f}\n")
    
    print("Final Welfare Equation:")
    print(f"{consumer_surplus:.4f} (Consumer Surplus) + {producer_surplus:.4f} (Producer Surplus) = {total_welfare:.4f} (Total Welfare)")
    
    # Final answer in the requested format
    print(f"\n<<<{total_welfare:.4f}>>>")

if __name__ == "__main__":
    main()