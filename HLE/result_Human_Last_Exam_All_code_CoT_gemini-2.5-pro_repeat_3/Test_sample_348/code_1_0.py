import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

def main():
    """
    This script calculates the market equilibrium and total welfare for the given
    supply and demand functions.
    """

    # 1. Define the supply and demand functions
    def demand_price(Q):
        """Demand function: P = 18 * e^(-arctan(Q))"""
        return 18 * np.exp(-np.arctan(Q))

    def supply_price(Q):
        """Supply function: P = ln(Q^3 - 2)"""
        # To handle potential domain errors during numerical solving,
        # return a large value if Q is outside the valid domain (Q^3 - 2 <= 0).
        with np.errstate(invalid='ignore'):
            if np.any(Q**3 - 2 <= 0):
                return np.inf
            return np.log(Q**3 - 2)

    # 2. Find the market equilibrium
    def equilibrium_equation(Q):
        """
        The equilibrium occurs where Demand Price = Supply Price.
        This function returns the difference, which we'll solve for zero.
        """
        return demand_price(Q) - supply_price(Q)

    # We need an initial guess to find the root. Let's try Q=4.
    initial_guess_Q = 4.0
    # Use fsolve to find the equilibrium quantity, Qe
    Q_e = fsolve(equilibrium_equation, initial_guess_Q)[0]

    # Calculate the equilibrium price, Pe
    P_e = demand_price(Q_e)

    # 3. Determine integration limits
    # The supply curve is defined for Q^3 - 2 > 0 => Q > 2^(1/3)
    lower_bound = 2**(1/3)

    # 4. Calculate Total Welfare
    # The integrand is the difference between the demand and supply curves
    integrand = equilibrium_equation
    
    # Use quad to perform the numerical integration
    total_welfare, _ = quad(integrand, lower_bound, Q_e)

    # Print the results step-by-step
    print("Given the functions:")
    print("Supply: P = ln(Q_S^3 - 2)")
    print("Demand: P = 18 * e^(-arctan(Q_D))\n")

    print("--- Step 1: Find Market Equilibrium ---")
    print("We solve for the quantity Q where Demand Price = Supply Price:")
    print("18 * e^(-arctan(Q)) = ln(Q^3 - 2)")
    print(f"Solving this numerically, we find the equilibrium quantity Qe = {Q_e:.4f}")
    print(f"The corresponding equilibrium price Pe = {P_e:.4f}\n")

    print("--- Step 2: Calculate Total Welfare ---")
    print("Total Welfare is the integral of (Demand Price - Supply Price).")
    print("The lower limit of integration is where the supply curve begins, Q = 2^(1/3).")
    print(f"Lower integration limit = 2^(1/3) = {lower_bound:.4f}")
    print("The final equation for total welfare is:")
    print(f"Total Welfare = ∫ (from {lower_bound:.4f} to {Q_e:.4f}) [ (18 * e^(-arctan(Q))) - (ln(Q^3 - 2)) ] dQ")
    print(f"\nThe calculated total welfare is: {total_welfare:.4f}")


if __name__ == "__main__":
    main()
<<<4.0722>>>