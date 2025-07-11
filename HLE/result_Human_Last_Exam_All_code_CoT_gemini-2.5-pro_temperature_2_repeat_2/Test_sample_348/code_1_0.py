import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

def main():
    """
    This script calculates the total welfare based on given supply and demand functions.
    """
    
    # 1. Define the supply and demand functions
    def supply_p(q):
        """Supply function: P = ln(Q^3 - 2)"""
        # Ensure the argument of log is positive to avoid math domain errors
        if q**3 - 2 <= 1e-9:
            return -np.inf 
        return np.log(q**3 - 2)

    def demand_p(q):
        """Demand function: P = 18 * e^(-arctan(Q))"""
        return 18 * np.exp(-np.arctan(q))

    # Function whose root is the equilibrium quantity
    def equilibrium_equation(q):
        q_val = q[0] # fsolve works with arrays
        # Ensure q is in the valid domain for supply
        if q_val**3 - 2 <= 0:
            # Return a large value if out of domain to guide the solver
            return 1e6
        return demand_p(q_val) - supply_p(q_val)

    # 2. Find Market Equilibrium
    print("Step 1: Finding the market equilibrium (Q_e, P_e).")
    
    # An initial guess for the root finding algorithm
    initial_guess = [4.5]
    
    # Solve for equilibrium quantity, Q_e
    q_e_solution = fsolve(equilibrium_equation, initial_guess)
    q_e = q_e_solution[0]
    
    # Calculate equilibrium price, P_e
    p_e = demand_p(q_e)
    
    print(f"Equilibrium Quantity (Q_e): {q_e:.4f}")
    print(f"Equilibrium Price (P_e): {p_e:.4f}\n")

    # 3. Determine Integration Bounds
    print("Step 2: Defining the integration bounds.")
    # For demand, the integral is from 0 to Q_e.
    # For supply, the function is defined for Q > 2^(1/3)
    q_min_supply = np.cbrt(2)
    print(f"Demand curve is integrated from 0 to Q_e.")
    print(f"Supply curve is defined for Q > 2^(1/3), so it is integrated from {q_min_supply:.4f} to Q_e.\n")

    # 4. Calculate Total Welfare
    print("Step 3: Calculating the components of total welfare.")
    
    # Integral of the demand curve from 0 to Q_e
    integral_demand, _ = quad(demand_p, 0, q_e)
    print(f"Area under Demand Curve (∫D(Q)dQ) from 0 to {q_e:.4f}: {integral_demand:.4f}")
    
    # Integral of the supply curve from Q_min to Q_e
    integral_supply, _ = quad(supply_p, q_min_supply, q_e)
    print(f"Area under Supply Curve (∫S(Q)dQ) from {q_min_supply:.4f} to {q_e:.4f}: {integral_supply:.4f}\n")

    # Total Welfare = Integral of Demand - Integral of Supply
    total_welfare = integral_demand - integral_supply
    
    # 5. Output the Final Result
    print("Step 4: Calculating final total welfare.")
    print("Total Welfare = ∫[0 to Q_e] D(Q) dQ - ∫[Q_min to Q_e] S(Q) dQ")
    print(f"Total Welfare = {integral_demand:.4f} - {integral_supply:.4f}")
    print(f"The total welfare is: {total_welfare:.4f}")

if __name__ == "__main__":
    main()
<<<11.4589>>>