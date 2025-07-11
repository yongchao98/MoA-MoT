import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

def main():
    """
    Calculates the total welfare given supply and demand functions.
    """
    
    # 1. Define the supply and demand functions
    def demand_p(Q):
        """Demand function P = 18 * e^(-arctan(Q))"""
        return 18 * np.exp(-np.arctan(Q))

    def supply_p(Q):
        """Supply function P = ln(Q^3 - 2)"""
        # The function is defined for Q^3 - 2 > 0, so Q > 2^(1/3)
        # We handle values outside the domain to assist the solver.
        with np.errstate(invalid='ignore', divide='ignore'):
            result = np.log(Q**3 - 2)
        # For Q outside the domain, log returns nan. Replace with a large value.
        if np.isscalar(result) and np.isnan(result):
             return np.inf
        return result

    # 2. Find the market equilibrium
    def equilibrium_equation(Q):
        """Equation to solve for equilibrium: Supply P - Demand P = 0"""
        return supply_p(Q) - demand_p(Q)

    # An initial guess for the solver, based on an approximate sketch of the curves.
    initial_guess_Q = 4.0
    # Use fsolve to find the equilibrium quantity Q*
    Q_star = fsolve(equilibrium_equation, initial_guess_Q)[0]
    # Calculate the equilibrium price P*
    P_star = demand_p(Q_star)

    # 3. Calculate Total Welfare
    
    # The lower bound of integration for the supply curve is where it is defined, i.e., Q^3-2 > 0
    Q_supply_min = (2)**(1/3)
    
    # Define the functions to be integrated (the same as the price functions)
    integrand_demand = demand_p
    integrand_supply = supply_p
    
    # Calculate the integral of the demand curve from 0 to Q*
    integral_demand, _ = quad(integrand_demand, 0, Q_star)
    
    # Calculate the integral of the supply curve from its start (Q_supply_min) to Q*
    integral_supply, _ = quad(integrand_supply, Q_supply_min, Q_star)
    
    # Total Welfare is the area under the demand curve minus the area under the supply curve
    total_welfare = integral_demand - integral_supply

    # 4. Print the results
    print("--- Market Equilibrium ---")
    print(f"Equilibrium Quantity (Q*): {Q_star}")
    print(f"Equilibrium Price (P*): {P_star}")
    print("\n--- Total Welfare Calculation ---")
    print(f"Total Welfare = Integral(Demand) - Integral(Supply)")
    print(f"Integral of Demand from 0 to Q*: {integral_demand}")
    print(f"Integral of Supply from {Q_supply_min:.4f} to Q*: {integral_supply}")
    print(f"\nFinal Total Welfare = {integral_demand:.4f} - {integral_supply:.4f}")
    print(f"Final Total Welfare = {total_welfare}")

if __name__ == "__main__":
    main()