import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

def solve_and_calculate_welfare():
    """
    This function calculates the market equilibrium and total welfare based on
    the given supply and demand functions.
    """
    # 1. Define supply and demand functions
    def supply(Q):
        # Supply function is valid only for Q^3 - 2 > 0.
        # This check prevents math domain errors during the solving process.
        if np.any(Q**3 <= 2):
            return np.inf  # Return a large number if outside the domain
        return np.log(Q**3 - 2)

    def demand(Q):
        return 18 * np.exp(-np.arctan(Q))

    # Function to find the root of (Demand - Supply = 0)
    def equilibrium_equation(Q):
        # fsolve works with arrays, so we extract the single element.
        q_val = Q[0]
        return demand(q_val) - supply(q_val)

    # 2. Find Market Equilibrium (Q_e, P_e)
    # Provide an initial guess for the solver. A simple plot or a few test values
    # suggest the equilibrium is around Q=4.5.
    initial_guess = [4.5]
    Q_e = fsolve(equilibrium_equation, initial_guess)[0]

    # Calculate equilibrium price using the demand function
    P_e = demand(Q_e)

    # 3. Calculate Consumer and Producer Surplus
    # The lower bound for supply integration is where the supply function is defined.
    # Q^3 - 2 > 0  => Q > 2^(1/3)
    min_q_supply = np.cbrt(2)
    
    # Integrate the demand function from 0 to Q_e
    integral_demand, _ = quad(demand, 0, Q_e)
    
    # Integrate the supply function from its minimum valid Q to Q_e
    integral_supply, _ = quad(supply, min_q_supply, Q_e)

    # Calculate consumer surplus
    consumer_surplus = integral_demand - (P_e * Q_e)
    
    # Calculate producer surplus
    producer_surplus = (P_e * Q_e) - integral_supply

    # 4. Calculate Total Welfare
    total_welfare = consumer_surplus + producer_surplus

    # Print the results as an equation
    print("Finding Total Welfare: Consumer Surplus + Producer Surplus")
    print(f"Equilibrium Quantity (Q_e) = {Q_e:.4f}")
    print(f"Equilibrium Price (P_e) = {P_e:.4f}")
    print("-" * 50)
    print(f"Consumer Surplus = {consumer_surplus:.4f}")
    print(f"Producer Surplus = {producer_surplus:.4f}")
    print("-" * 50)
    print("Total Welfare = Consumer Surplus + Producer Surplus")
    print(f"Total Welfare = {consumer_surplus:.4f} + {producer_surplus:.4f}")
    print(f"Total Welfare = {total_welfare:.4f}")

# Run the calculation and print the output
solve_and_calculate_welfare()