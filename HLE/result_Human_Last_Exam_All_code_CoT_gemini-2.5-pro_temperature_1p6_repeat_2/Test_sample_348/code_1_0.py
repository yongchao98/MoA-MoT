import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

def solve_total_welfare():
    """
    This script calculates the total welfare by finding the market equilibrium,
    computing consumer and producer surplus, and summing them.
    """

    # --- 1. Define Supply and Demand Functions ---

    def demand(Q):
        """Demand function P = 18 * e^(-arctan(Q))"""
        return 18 * np.exp(-np.arctan(Q))

    def supply(Q):
        """Supply function P = ln(Q^3 - 2)"""
        # np.log will return nan for non-positive inputs, which is handled.
        return np.log(Q**3 - 2)

    # --- 2. Find Market Equilibrium ---

    def equilibrium_equation(Q):
        """Equation to solve for equilibrium: D(Q) - S(Q) = 0"""
        q_val = Q[0]
        # Supply is undefined for q_val <= 2^(1/3). Return a large number
        # to guide the solver away from this invalid region.
        if q_val**3 <= 2:
            return 1e6
        return demand(q_val) - supply(q_val)

    # Solve for equilibrium quantity Q*
    initial_guess = [4.5]
    Q_equilibrium = fsolve(equilibrium_equation, initial_guess)[0]

    # Calculate equilibrium price P*
    P_equilibrium = demand(Q_equilibrium)

    # --- 3. Calculate Consumer and Producer Surplus ---

    # Consumer Surplus (CS) = ∫[0 to Q*] (D(Q) - P*) dQ
    integrand_cs = lambda Q: demand(Q) - P_equilibrium
    consumer_surplus, _ = quad(integrand_cs, 0, Q_equilibrium)

    # Producer Surplus (PS) = ∫[Q_min to Q*] (P* - S(Q)) dQ
    # The supply curve is defined for Q > 2^(1/3)
    Q_min_supply = (2.0)**(1/3.0)
    integrand_ps = lambda Q: P_equilibrium - supply(Q)
    producer_surplus, _ = quad(integrand_ps, Q_min_supply, Q_equilibrium)

    # --- 4. Calculate Total Welfare and Print Results ---

    total_welfare = consumer_surplus + producer_surplus

    print("Step 1: Find Market Equilibrium")
    print(f"The equilibrium quantity (Q*) is found by solving ln(Q^3 - 2) = 18 * e^(-arctan(Q)).")
    print(f"Equilibrium Quantity (Q*) = {Q_equilibrium:.4f}")
    print(f"Equilibrium Price (P*) = {P_equilibrium:.4f}\n")

    print("Step 2: Calculate Consumer Surplus (CS)")
    print(f"CS = integral from 0 to {Q_equilibrium:.4f} of [18 * e^(-arctan(Q)) - {P_equilibrium:.4f}] dQ")
    print(f"Consumer Surplus = {consumer_surplus:.4f}\n")

    print("Step 3: Calculate Producer Surplus (PS)")
    print(f"The supply curve starts where Q > 2^(1/3), so Q_min = {Q_min_supply:.4f}")
    print(f"PS = integral from {Q_min_supply:.4f} to {Q_equilibrium:.4f} of [{P_equilibrium:.4f} - ln(Q^3 - 2)] dQ")
    print(f"Producer Surplus = {producer_surplus:.4f}\n")

    print("Step 4: Calculate Total Welfare (TW)")
    print("Total Welfare = Consumer Surplus + Producer Surplus")
    print(f"Total Welfare = {consumer_surplus:.4f} + {producer_surplus:.4f}")
    print(f"Total Welfare = {total_welfare:.4f}")
    
    # Final answer for automated evaluation
    # print(f"\n<<<{total_welfare:.4f}>>>")

solve_total_welfare()