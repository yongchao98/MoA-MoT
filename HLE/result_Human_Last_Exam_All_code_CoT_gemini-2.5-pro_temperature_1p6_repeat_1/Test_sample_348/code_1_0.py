import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad

def main():
    """
    Calculates the market equilibrium and total welfare for the given supply and demand functions.
    """

    # Define the supply and demand functions
    def demand_price(q):
        """Calculates price from the demand function for a given quantity q."""
        return 18 * np.exp(-np.arctan(q))

    def supply_price(q):
        """Calculates price from the supply function for a given quantity q."""
        # The supply function is only defined for q^3 - 2 > 0.
        # This handles arrays and single values to avoid math domain errors.
        with np.errstate(all='ignore'):
            return np.log(q**3 - 2)

    # Define the equation to find equilibrium: Demand Price - Supply Price = 0
    def equilibrium_equation(q):
        """Difference between demand and supply price. The root is the equilibrium quantity."""
        # fsolve passes an array, so we access the first element.
        q_val = q[0]
        # Return a large number if outside the supply domain to guide the solver.
        if q_val**3 - 2 <= 0:
            return 1e6
        return demand_price(q_val) - supply_price(q_val)

    # --- Step 1: Find Market Equilibrium ---
    
    # An initial guess for the solver. It must be greater than ∛2 ≈ 1.26.
    initial_guess = [2.0]
    
    # Use fsolve to find the equilibrium quantity (Q_E)
    q_equilibrium_array = fsolve(equilibrium_equation, initial_guess)
    q_equilibrium = q_equilibrium_array[0]
    
    # Calculate the equilibrium price (P_E) using the demand function
    p_equilibrium = demand_price(q_equilibrium)

    # --- Step 2: Calculate Total Welfare ---

    # The lower bound of integration for the supply function is where it is first defined.
    q_min_supply = (2)**(1/3)
    
    # Calculate the integral of the demand curve from 0 to Q_E
    integral_demand, _ = quad(demand_price, 0, q_equilibrium)
    
    # Calculate the integral of the supply curve from its defined start to Q_E
    integral_supply, _ = quad(supply_price, q_min_supply, q_equilibrium)
    
    # Calculate the final total welfare
    total_welfare = integral_demand - integral_supply

    # --- Step 3: Output the Results ---
    
    print("Market Equilibrium Calculation:")
    print(f"Equilibrium Quantity (Q_E): {q_equilibrium:.4f}")
    print(f"Equilibrium Price (P_E): {p_equilibrium:.4f}")
    print("-" * 40)
    
    print("Total Welfare Calculation:")
    print("Total Welfare = ∫Demand(Q)dQ - ∫Supply(Q)dQ")
    print(f"Demand Integral = ∫ from 0 to {q_equilibrium:.4f} (18e^(-arctan(Q))) dQ")
    print(f"Supply Integral = ∫ from {q_min_supply:.4f} to {q_equilibrium:.4f} (ln(Q³-2)) dQ")
    print("-" * 40)

    # Output the final calculation with each number
    print("Final Equation:")
    print(f"{integral_demand:.4f} - {integral_supply:.4f} = {total_welfare:.4f}")
    
    # Return the final numerical answer for extraction
    return total_welfare

if __name__ == "__main__":
    main()
