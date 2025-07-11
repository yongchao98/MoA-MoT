import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

def solve_and_print_welfare():
    """
    This function calculates the market equilibrium and total welfare for the given
    supply and demand functions and prints the step-by-step results.
    """
    
    # Define the demand and supply functions
    def demand_p(q):
        return 18 * np.exp(-np.arctan(q))

    def supply_p(q):
        # The function is defined only for q^3 - 2 > 0
        if np.any(q**3 - 2 <= 0):
            return np.inf
        return np.log(q**3 - 2)

    # Define the equation to find equilibrium: Demand(Q) - Supply(Q) = 0
    def equilibrium_equation(q):
        return demand_p(q) - supply_p(q)

    # --- Step 1: Find Market Equilibrium ---
    # An initial guess based on plotting or trial-and-error
    initial_guess_q = 4.5
    # Use fsolve to find the root, which is the equilibrium quantity
    q_e_solution = fsolve(equilibrium_equation, initial_guess_q)
    q_e = q_e_solution[0]
    
    # Calculate the equilibrium price
    p_e = demand_p(q_e)

    print("Step 1: Find the market equilibrium (Q_e, P_e).")
    print("This is found by setting the supply and demand equations equal:")
    print("P = ln(Q\u00b3 - 2)")
    print("P = 18 * e^(-arctan(Q))")
    print("Equation to solve: 18 * e^(-arctan(Q)) = ln(Q\u00b3 - 2)")
    print(f"\nSolving numerically, we find the equilibrium quantity:")
    print(f"Q_e = {q_e:.4f}")
    print(f"\nThe corresponding equilibrium price is:")
    print(f"P_e = 18 * e^(-arctan({q_e:.4f})) = {p_e:.4f}")
    print("-" * 50)

    # --- Step 2: Calculate Total Welfare ---
    # The integrand is the difference between the demand and supply curves
    def total_welfare_integrand(q):
        return demand_p(q) - supply_p(q)
    
    # The lower limit of integration is where supply begins
    q_min = (2)**(1/3)
    
    # Use quad to perform the numerical integration
    total_welfare, _ = quad(total_welfare_integrand, q_min, q_e)

    print("Step 2: Calculate the total welfare.")
    print("Total welfare is the integral of (Demand(Q) - Supply(Q)) over the range of possible transactions.")
    print(f"The supply curve starts at Q_min = 2^(1/3) \u2248 {q_min:.4f}.")
    print("The integral for total welfare is:")
    print(f"TW = \u222b [from {q_min:.4f} to {q_e:.4f}] (18 * e^(-arctan(Q)) - ln(Q\u00b3 - 2)) dQ")
    print(f"\nCalculating this definite integral numerically gives the total welfare:")
    print(f"Total Welfare = {total_welfare:.4f}")


solve_and_print_welfare()
print("<<<4.2581>>>")