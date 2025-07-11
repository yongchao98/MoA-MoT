import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve

def solve_and_print_welfare():
    """
    This function calculates and prints the total welfare based on the given
    supply and demand functions.
    """
    # Define the demand function P(Q)
    def demand(q):
        return 18 * np.exp(-np.arctan(q))

    # Define the supply function P(Q)
    # The function is defined for Q^3 - 2 > 0.
    def supply(q):
        # The solver might test values where the log is undefined.
        # We return a large number to push the solver away from this region.
        with np.errstate(invalid='ignore'): # Suppress log(negative) warnings
            val = np.log(q**3 - 2)
        if np.isnan(val):
            return 1e9 # Return a large number for invalid inputs
        return val

    # Define the equation to find the equilibrium quantity: Demand(Q) - Supply(Q) = 0
    def equilibrium_equation(q):
        return demand(q) - supply(q)

    # --- Step 1: Find the Market Equilibrium Quantity (Q_E) ---
    # An initial guess is needed for the numerical solver.
    # A quick plot or trial-and-error suggests the intersection is between 4 and 5.
    initial_guess = 4.5
    # Use fsolve to find the root of the equilibrium equation.
    Q_E = fsolve(equilibrium_equation, initial_guess)[0]

    # --- Step 2: Define Integration Bounds ---
    # The supply curve is defined for Q^3 - 2 > 0, so Q > 2^(1/3)
    supply_q_min = 2**(1/3)

    # --- Step 3: Calculate the Integrals for Total Welfare ---
    # Total Welfare = ∫[0 to Q_E] D(Q)dQ - ∫[Q_min to Q_E] S(Q)dQ
    
    # Integral of the demand curve from 0 to Q_E
    integral_demand, _ = quad(demand, 0, Q_E)

    # Integral of the supply curve from its starting point to Q_E
    integral_supply, _ = quad(supply, supply_q_min, Q_E)

    # --- Step 4: Calculate Total Welfare and Print Results ---
    total_welfare = integral_demand - integral_supply
    
    print("To find the total welfare, we use the formula:")
    print("Total Welfare = ∫[0 to Q_E] D(Q)dQ - ∫[Q_min to Q_E] S(Q)dQ\n")
    
    print(f"1. The equilibrium quantity (Q_E) where demand equals supply is: {Q_E:.4f}")
    print(f"2. The minimum quantity for the supply curve (Q_min) is 2^(1/3): {supply_q_min:.4f}\n")
    
    print("3. Calculating the definite integrals:")
    print(f"   - Integral of Demand from 0 to {Q_E:.4f} = {integral_demand:.4f}")
    print(f"   - Integral of Supply from {supply_q_min:.4f} to {Q_E:.4f} = {integral_supply:.4f}\n")
    
    print("4. Final Calculation for Total Welfare:")
    print(f"   Total Welfare = {integral_demand:.4f} - {integral_supply:.4f}")
    print(f"   Total Welfare = {total_welfare:.4f}")

# Run the calculation and print the output
solve_and_print_welfare()
<<<18.3308>>>