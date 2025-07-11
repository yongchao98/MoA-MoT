import numpy as np
from scipy.optimize import root_scalar
from scipy.integrate import quad

def main():
    """
    Calculates and prints the total welfare based on given supply and demand functions.
    """
    # Define the supply and demand functions
    # Supply: P = ln(Q^3 - 2)
    def supply_func(q):
        # The argument of ln must be > 0. Q^3 - 2 > 0 => Q > 2^(1/3)
        # We add a small epsilon to avoid issues at the boundary for the solver.
        if np.any(q**3 - 2 <= 1e-9):
            return np.inf # Return a large value if Q is out of domain
        return np.log(q**3 - 2)

    # Demand: P = 18 * e^(-arctan(Q))
    def demand_func(q):
        return 18 * np.exp(-np.arctan(q))

    # --- 1. Find Market Equilibrium ---
    # We need to find the root of Demand(Q) - Supply(Q) = 0
    def equilibrium_equation(q):
        return demand_func(q) - supply_func(q)

    # By plotting or testing values, we can see the root is between Q=4 and Q=5.
    try:
        solution = root_scalar(equilibrium_equation, bracket=[4, 5])
        if not solution.converged:
            print("Numerical solver failed to converge.")
            return
        Q_e = solution.root
    except ValueError as e:
        print(f"Error finding equilibrium quantity: {e}")
        return

    # --- 2. Calculate Equilibrium Price ---
    P_e = demand_func(Q_e)

    # --- 3. Calculate Consumer Surplus (CS) ---
    # CS = Integral(Demand(Q) dQ) from 0 to Q_e - P_e * Q_e
    integral_demand, _ = quad(demand_func, 0, Q_e)
    consumer_surplus = integral_demand - (P_e * Q_e)

    # --- 4. Calculate Producer Surplus (PS) ---
    # PS = P_e * Q_e - Integral(Supply(Q) dQ) from Q_min to Q_e
    # Find Q_min: Q^3 - 2 > 0 => Q > 2^(1/3)
    Q_min = 2**(1/3)
    integral_supply, _ = quad(supply_func, Q_min, Q_e)
    producer_surplus = (P_e * Q_e) - integral_supply

    # --- 5. Calculate Total Welfare ---
    total_welfare = consumer_surplus + producer_surplus

    # Print the final equation with the calculated values
    print("Final Calculation:")
    print(f"Total Welfare = Consumer Surplus + Producer Surplus")
    print(f"Total Welfare = {consumer_surplus:.4f} + {producer_surplus:.4f} = {total_welfare:.4f}")

if __name__ == "__main__":
    main()