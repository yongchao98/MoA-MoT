import numpy as np
from scipy.optimize import root_scalar
from scipy.integrate import quad

def solve_and_print_welfare():
    """
    This function calculates and prints the market equilibrium and total welfare
    for the given supply and demand functions.
    """
    # Step 1: Define the supply and demand functions
    def demand(q):
        """Demand function: P = 18 * e^(-arctan(Q))"""
        return 18 * np.exp(-np.arctan(q))

    def supply(q):
        """Supply function: P = ln(Q^3 - 2)"""
        # This function is only defined for q > 2^(1/3)
        # The solver will keep q in the valid domain using a bracket.
        return np.log(q**3 - 2)

    # Define the equation whose root is the equilibrium quantity
    def equilibrium_equation(q):
        """Returns Demand(Q) - Supply(Q). Root is at equilibrium."""
        return demand(q) - supply(q)

    # Step 2: Find the market equilibrium
    # The supply function is defined for Q^3 - 2 > 0, so Q > 2^(1/3)
    q_min = 2**(1/3)

    # We use a numerical solver to find the equilibrium quantity Q_e
    # by finding the root of equilibrium_equation(Q) = 0.
    # We search in a bracket safely inside the supply function's domain.
    # A quick analysis shows the root is likely between 2 and 10.
    try:
        sol = root_scalar(equilibrium_equation, bracket=[q_min + 1e-9, 10])
        q_equilibrium = sol.root
        # Calculate the equilibrium price P_e
        p_equilibrium = demand(q_equilibrium)
    except ValueError as e:
        print(f"An error occurred during equilibrium calculation: {e}")
        print("Could not find an equilibrium point in the given bracket.")
        return

    # Step 3: Calculate the Total Welfare using the formula:
    # TW = Integral[0 to Q_e](Demand(Q))dQ - Integral[Q_min to Q_e](Supply(Q))dQ - P_e * Q_min

    # Calculate the integral of the demand function from 0 to Q_e
    integral_demand, _ = quad(demand, 0, q_equilibrium)

    # Calculate the integral of the supply function from Q_min to Q_e
    integral_supply, _ = quad(supply, q_min, q_equilibrium)

    # Calculate the extra term from the producer surplus formula
    pe_qmin_term = p_equilibrium * q_min

    # Calculate the final total welfare
    total_welfare = integral_demand - integral_supply - pe_qmin_term

    # Step 4: Output the results clearly
    print("--- Market Equilibrium Calculation ---")
    print(f"Supply function: P = ln(Q^3 - 2)")
    print(f"Demand function: P = 18 * e^(-arctan(Q))")
    print(f"Equilibrium is found where Demand(Q) = Supply(Q).")
    print(f"Solved numerically, the equilibrium quantity Q_e is: {q_equilibrium:.4f}")
    print(f"The corresponding equilibrium price P_e is: {p_equilibrium:.4f}")
    
    print("\n--- Total Welfare Calculation ---")
    print("Total Welfare is the area between the demand and supply curves.")
    print("The formula for Total Welfare in this case is:")
    print("TW = Integral[0 to Q_e](Demand(Q))dQ - (Integral[Q_min to Q_e](Supply(Q))dQ + P_e * Q_min)")
    
    print(f"\nWhere:")
    print(f"Q_e = {q_equilibrium:.4f}")
    print(f"P_e = {p_equilibrium:.4f}")
    print(f"Q_min (supply curve start) = 2^(1/3) = {q_min:.4f}")
    
    print("\nCalculating the components:")
    print(f"1. Integral of Demand from 0 to {q_equilibrium:.4f} = {integral_demand:.4f}")
    print(f"2. Integral of Supply from {q_min:.4f} to {q_equilibrium:.4f} = {integral_supply:.4f}")
    print(f"3. P_e * Q_min term = {p_equilibrium:.4f} * {q_min:.4f} = {pe_qmin_term:.4f}")
    
    print("\n--- Final Equation with Values ---")
    print(f"Total Welfare = {integral_demand:.4f} - {integral_supply:.4f} - {pe_qmin_term:.4f}")
    print(f"Total Welfare = {total_welfare:.4f}")

# Execute the function to get the answer
solve_and_print_welfare()