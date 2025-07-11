import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

def main():
    """
    This script calculates the total welfare based on given supply and demand functions.
    """
    
    # Step 1: Define the supply and demand functions
    def supply(q):
        """Supply function: P = ln(Q^3 - 2)"""
        # Function is defined for Q^3 - 2 > 0.
        # Handle cases outside the domain to avoid math errors.
        if q**3 <= 2:
            return -np.inf  # Use a value that shows it's outside the valid economic range
        return np.log(q**3 - 2)

    def demand(q):
        """Demand function: P = 18 * e^(-arctan(Q))"""
        return 18 * np.exp(-np.arctan(q))

    # Step 2: Find the market equilibrium
    # We need to solve Demand(Q) = Supply(Q)
    def equilibrium_equation(q):
        # We need to ensure q is in the domain of the supply function for the solver.
        # fsolve can sometimes try values outside the domain, causing errors.
        if q**3 <= 2:
            # Return a large value to push the solver towards the valid domain
            return 1e9
        return demand(q) - supply(q)

    # Use an initial guess for the solver. Based on manual checking, Q is between 4 and 5.
    initial_guess = 4.0
    # Use fsolve to find the root, which is the equilibrium quantity (Q_e)
    equilibrium_quantity = fsolve(equilibrium_equation, initial_guess)[0]
    
    # Calculate the equilibrium price (P_e) using the demand function
    equilibrium_price = demand(equilibrium_quantity)
    
    # Step 3: Determine integration bounds and calculate Total Welfare
    # Total Welfare = Integral(Demand(Q))dQ - Integral(Supply(Q))dQ
    
    # The lower bound for supply is where P >= 0, so ln(Q^3 - 2) >= 0 => Q^3 >= 3
    supply_lower_bound = (3)**(1/3)
    
    # The lower bound for demand is 0
    demand_lower_bound = 0
    
    # Calculate the integral for the area under the demand curve from 0 to Q_e
    integral_demand, _ = quad(demand, demand_lower_bound, equilibrium_quantity)
    
    # Calculate the integral for the area under the supply curve from its start to Q_e
    integral_supply, _ = quad(supply, supply_lower_bound, equilibrium_quantity)
    
    # Total welfare is the difference
    total_welfare = integral_demand - integral_supply
    
    # Step 4: Print the results step-by-step
    print("Finding Market Equilibrium...")
    print(f"Equilibrium Quantity (Q_e) = {equilibrium_quantity:.4f}")
    print(f"Equilibrium Price (P_e) = {equilibrium_price:.4f}\n")
    
    print("Calculating Total Welfare...")
    print("Formula: Total Welfare = Integral(Demand) - Integral(Supply)")
    print(f"Integral of Demand from {demand_lower_bound} to {equilibrium_quantity:.4f} = {integral_demand:.4f}")
    print(f"Integral of Supply from {supply_lower_bound:.4f} to {equilibrium_quantity:.4f} = {integral_supply:.4f}\n")
    
    print("Final Calculation:")
    print(f"Total Welfare = {integral_demand:.4f} - {integral_supply:.4f}")
    print(f"Total Welfare = {total_welfare:.4f}")

if __name__ == "__main__":
    main()
<<<18.1565>>>