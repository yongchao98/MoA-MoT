import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

def main():
    """
    This script calculates the total welfare based on given supply and demand functions.
    """
    # 1. Define Supply and Demand Functions
    # Supply: P = ln(Q^3 - 2)
    # This function is defined for Q > 2^(1/3)
    def supply_price(q):
        # Adding a safeguard for the domain, though quad integration bounds should prevent issues.
        if q**3 - 2 <= 1e-9: # Use a small epsilon for floating point safety
            return -np.inf # Should not be reached with correct integration limits
        return np.log(q**3 - 2)

    # Demand: P = 18 * e^(-arctan(Q))
    def demand_price(q):
        return 18 * np.exp(-np.arctan(q))

    # 2. Find Market Equilibrium
    # We need to solve P_demand(Q) = P_supply(Q)
    # Define the equation for the solver: demand_price(Q) - supply_price(Q) = 0
    # fsolve expects a function that takes a list/array as input
    def equilibrium_equation(q_arr):
        q = q_arr[0]
        return demand_price(q) - supply_price(q)

    # The supply curve is defined for Q > 2^(1/3) ~= 1.26
    # We must provide an initial guess within this domain.
    initial_guess_q = [2.0]
    
    # Use fsolve to find the equilibrium quantity Q_E
    q_equilibrium_solution = fsolve(equilibrium_equation, initial_guess_q)
    q_equilibrium = q_equilibrium_solution[0]

    # 3. Calculate Total Welfare
    # Total Welfare = Integral(Demand) - Integral(Supply) over the relevant ranges

    # The integral under the demand curve from 0 to Q_E
    integral_demand, _ = quad(demand_price, 0, q_equilibrium)

    # The lower bound for the supply integral is where the supply curve begins.
    # P = ln(Q^3 - 2) is defined for Q^3 - 2 > 0 => Q > 2^(1/3)
    q_min_supply = 2**(1/3)

    # The integral under the supply curve from its start (q_min_supply) to Q_E
    integral_supply, _ = quad(supply_price, q_min_supply, q_equilibrium)

    # Total welfare is the difference between the two integrals.
    total_welfare = integral_demand - integral_supply

    # 4. Print the final results, showing each number in the equation.
    print("The final calculation for Total Welfare is based on the formula:")
    print("Total Welfare = ∫[0 to Q_E] P_Demand(Q) dQ - ∫[Q_min to Q_E] P_Supply(Q) dQ\n")
    print(f"Equilibrium Quantity (Q_E) found at: {q_equilibrium:.4f}")
    print(f"Integral of Demand from 0 to {q_equilibrium:.4f} = {integral_demand:.4f}")
    print(f"Integral of Supply from {q_min_supply:.4f} to {q_equilibrium:.4f} = {integral_supply:.4f}\n")
    
    print("Final Equation:")
    print(f"Total Welfare = {integral_demand:.4f} - {integral_supply:.4f} = {total_welfare:.4f}")

if __name__ == "__main__":
    main()