import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

def find_total_welfare():
    """
    This function calculates the market equilibrium and total welfare
    for the given supply and demand functions.
    """
    # 1. Define the supply and demand price functions
    def supply_price(Q):
        # The argument of ln must be positive, otherwise, it's undefined.
        # We handle this by returning a very large number to guide the solver.
        if Q**3 - 2 <= 0:
            return np.inf
        return np.log(Q**3 - 2)

    def demand_price(Q):
        return 18 * np.exp(-np.arctan(Q))

    # 2. Find the market equilibrium (Q_e, P_e)
    # We define a function whose root is the equilibrium quantity.
    # Equilibrium occurs where supply_price(Q) - demand_price(Q) = 0.
    def equilibrium_equation(Q):
        return supply_price(Q) - demand_price(Q)

    # We need an initial guess for the root-finding algorithm.
    # The supply function requires Q^3 > 2, so Q > 2^(1/3) ≈ 1.26.
    # Let's try an initial guess slightly above this, e.g., 4.0.
    initial_guess_Q = 4.0
    # Use fsolve to find the equilibrium quantity, Q_e. fsolve returns an array.
    Q_e_solution = fsolve(equilibrium_equation, initial_guess_Q)
    Q_e = Q_e_solution[0]

    # Calculate the equilibrium price, P_e, using the demand function.
    P_e = demand_price(Q_e)

    # 3. Calculate Total Welfare
    # The integrand is the difference between the demand and supply price curves.
    def welfare_integrand(Q):
        return demand_price(Q) - supply_price(Q)

    # The lower bound of integration is the minimum quantity where supply is defined.
    lower_integration_bound = (2)**(1/3)

    # Use quad to calculate the definite integral for total welfare.
    # The upper bound is the equilibrium quantity, Q_e.
    total_welfare, integration_error = quad(welfare_integrand, lower_integration_bound, Q_e)

    # 4. Print the results
    print("Step 1: Find Market Equilibrium")
    print(f"Equilibrium Quantity (Q_e) = {Q_e:.4f}")
    print(f"Equilibrium Price (P_e) = {P_e:.4f}\n")

    print("Step 2: Calculate Total Welfare")
    print("Total Welfare is the integral of (Demand Price - Supply Price)")
    print(f"The equation to solve is:")
    print(f"Integral from Q_min={lower_integration_bound:.4f} to Q_e={Q_e:.4f} of [18*exp(-arctan(Q)) - ln(Q^3 - 2)] dQ")
    print(f"\nCalculated Total Welfare = {total_welfare:.4f}")

    # Return the final numerical answer for the platform
    return total_welfare

if __name__ == '__main__':
    final_answer = find_total_welfare()
    # print(f"<<<{final_answer}>>>") # Final answer tag for the system

# To directly run and see the output
find_total_welfare()