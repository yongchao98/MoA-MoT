import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

def solve_and_print_welfare():
    """
    This function calculates and prints the total welfare based on given
    supply and demand functions.
    """
    # 1. Define the supply and demand functions
    def demand_func(q):
        """Demand function: P = 18 * e^(-arctan(Q))"""
        return 18 * np.exp(-np.arctan(q))

    def supply_func(q):
        """Supply function: P = ln(Q^3 - 2)"""
        # The function is defined for Q^3 - 2 > 0.
        # Handle values outside the domain by returning a large number to aid the solver.
        if q**3 <= 2:
            return np.inf
        return np.log(q**3 - 2)

    # 2. Find market equilibrium by solving Demand(Q) = Supply(Q)
    def equilibrium_equation(q):
        """The equation to solve for equilibrium: Demand(Q) - Supply(Q) = 0"""
        return demand_func(q) - supply_func(q)

    # Initial guess for the solver, based on plotting or manual checks
    initial_guess = 4.5
    # Find the equilibrium quantity (Qe)
    q_equilibrium = fsolve(equilibrium_equation, initial_guess)[0]

    # Calculate the equilibrium price (Pe)
    p_equilibrium = demand_func(q_equilibrium)
    
    # 3. Calculate Total Welfare
    # The minimum quantity for the supply function to be defined is Q > 2^(1/3)
    q_min_supply = np.cbrt(2)

    # Integrate the demand function from 0 to Qe
    integral_demand, _ = quad(demand_func, 0, q_equilibrium)

    # Integrate the supply function from its minimum defined Q to Qe
    integral_supply, _ = quad(supply_func, q_min_supply, q_equilibrium)
    
    # Total Welfare is the area between the curves
    total_welfare = integral_demand - integral_supply

    # --- Outputting the results ---
    print("### Market Welfare Calculation ###\n")
    
    print("1. Find Market Equilibrium (Qe, Pe):")
    print("   We solve for Q where Supply equals Demand:")
    print(f"   ln(Q^3 - 2) = 18 * e^(-arctan(Q))\n")
    print(f"   Using a numerical solver, the equilibrium quantity is:")
    print(f"   Qe = {q_equilibrium:.4f}\n")
    print("   The equilibrium price Pe is found by substituting Qe into the demand function:")
    print(f"   Pe = 18 * e^(-arctan({q_equilibrium:.4f}))")
    print(f"   Pe = {p_equilibrium:.4f}\n")

    print("2. Calculate Total Welfare:")
    print("   Total Welfare = (Integral of Demand) - (Integral of Supply)\n")
    
    print("   Integral of Demand from 0 to Qe:")
    print(f"   ∫[0 to {q_equilibrium:.4f}] (18 * e^(-arctan(Q))) dQ = {integral_demand:.4f}\n")

    print("   Integral of Supply from its start point (2^(1/3)) to Qe:")
    print(f"   (Supply is defined for Q > {q_min_supply:.4f})")
    print(f"   ∫[{q_min_supply:.4f} to {q_equilibrium:.4f}] (ln(Q^3 - 2)) dQ = {integral_supply:.4f}\n")

    print("   Final Calculation:")
    print(f"   Total Welfare = {integral_demand:.4f} - {integral_supply:.4f} = {total_welfare:.4f}")

    print(f"\n<<<{total_welfare:.4f}>>>")

solve_and_print_welfare()