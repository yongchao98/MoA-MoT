import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

def main():
    """
    This script calculates the total welfare based on given supply and demand functions.
    """
    # Define the demand function
    def demand(Q_D):
        """Demand function: P = 18 * e^(-arctan(Q_D))"""
        return 18 * np.exp(-np.arctan(Q_D))

    # Define the supply function
    def supply(Q_S):
        """Supply function: P = ln(Q_S^3 - 2)"""
        # The function is defined for Q_S^3 - 2 > 0.
        # We add a condition to handle values outside the domain for the solver.
        if np.any(Q_S**3 <= 2):
            return np.inf
        return np.log(Q_S**3 - 2)

    # --- Step 1: Find the Market Equilibrium ---
    print("Step 1: Finding the market equilibrium point (Q_e, P_e).")
    
    # Define the equation to find the equilibrium quantity where Demand(Q) = Supply(Q)
    def equilibrium_equation(Q):
        return demand(Q) - supply(Q)

    # Use a numerical solver (fsolve) to find the root, which is Q_e.
    # An initial guess of 4.5 is used based on preliminary analysis.
    initial_guess_Q = 4.5
    Q_e_solution = fsolve(equilibrium_equation, x0=initial_guess_Q)
    Q_e = Q_e_solution[0]

    # Calculate the equilibrium price (P_e) using the demand function
    P_e = demand(Q_e)
    
    print(f"   Equilibrium Quantity (Q_e): {Q_e:.4f}")
    print(f"   Equilibrium Price (P_e): {P_e:.4f}\n")

    # --- Step 2: Calculate Consumer Surplus (CS) ---
    print("Step 2: Calculating Consumer Surplus (CS).")
    # CS = Integral of Demand(Q) from 0 to Q_e - (P_e * Q_e)
    integral_demand, _ = quad(demand, 0, Q_e)
    consumer_surplus = integral_demand - (P_e * Q_e)
    print(f"   Consumer Surplus = (Integral of Demand from 0 to {Q_e:.2f}) - ({P_e:.2f} * {Q_e:.2f})")
    print(f"   Consumer Surplus = {integral_demand:.4f} - {P_e * Q_e:.4f} = {consumer_surplus:.4f}\n")

    # --- Step 3: Calculate Producer Surplus (PS) ---
    print("Step 3: Calculating Producer Surplus (PS).")
    # The supply curve starts where Q^3 - 2 > 0, so Q_min = 2^(1/3)
    Q_min = np.cbrt(2)
    # PS = (P_e * Q_e) - Integral of Supply(Q) from Q_min to Q_e
    integral_supply, _ = quad(supply, Q_min, Q_e)
    producer_surplus = (P_e * Q_e) - integral_supply
    print(f"   Producer Surplus = ({P_e:.2f} * {Q_e:.2f}) - (Integral of Supply from {Q_min:.2f} to {Q_e:.2f})")
    print(f"   Producer Surplus = {P_e * Q_e:.4f} - {integral_supply:.4f} = {producer_surplus:.4f}\n")

    # --- Step 4: Calculate Total Welfare (TW) ---
    print("Step 4: Calculating Total Welfare (TW).")
    total_welfare = consumer_surplus + producer_surplus
    
    # Final output showing the equation
    print("Final Calculation:")
    print(f"Total Welfare = Consumer Surplus + Producer Surplus")
    print(f"Total Welfare = {consumer_surplus:.4f} + {producer_surplus:.4f} = {total_welfare:.4f}")

if __name__ == "__main__":
    main()