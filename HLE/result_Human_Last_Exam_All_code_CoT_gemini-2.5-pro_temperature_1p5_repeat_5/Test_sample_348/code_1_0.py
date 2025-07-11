import numpy as np
from scipy.optimize import root_scalar
from scipy.integrate import quad

def main():
    """
    Calculates the total welfare by finding market equilibrium,
    consumer surplus, and producer surplus.
    """
    # 1. Define supply and demand functions
    def demand_p(q):
        """Demand function: P = 18 * e^(-arctan(Q))"""
        # Ensure q is non-negative
        q = np.maximum(q, 0)
        return 18 * np.exp(-np.arctan(q))

    def supply_p(q):
        """Supply function: P = ln(Q^3 - 2)"""
        # The function is defined for Q^3 - 2 > 0, which means Q > ∛2.
        # Handle the domain for the solver and integrator.
        if isinstance(q, (float, int)):
            if q**3 <= 2:
                # Return a value that indicates an invalid domain,
                # e.g., a large number for a minimization problem
                # or handle appropriately in the context of the solver.
                # Here we use it for a root finder where P_demand - P_supply = 0.
                # Returning a very low price (negative infinity) ensures
                # demand > supply well below equilibrium.
                return -np.inf
        
        # For array inputs, handle element-wise
        q_safe = np.maximum(q, np.cbrt(2) + 1e-9) # Add epsilon to avoid log(0)
        return np.log(q_safe**3 - 2)

    def equilibrium_equation(q):
        """Difference between demand and supply, for finding the root."""
        return demand_p(q) - supply_p(q)

    # 2. Find market equilibrium
    # The supply curve is defined for Q > ∛2 ≈ 1.26.
    # By plotting or testing values, we can see the intersection is between Q=4 and Q=5.
    # equilibrium_equation(4) ≈ 0.65 > 0
    # equilibrium_equation(5) ≈ -0.26 < 0
    # So a root exists in the bracket [4, 5].
    q_min_supply = np.cbrt(2)
    try:
        solution = root_scalar(equilibrium_equation, bracket=[q_min_supply + 0.1, 10])
        Qe = solution.root
        Pe = demand_p(Qe)
    except ValueError as e:
        print(f"Error finding equilibrium: {e}")
        return

    # 3. Calculate Consumer Surplus (CS)
    # CS = Integral from 0 to Qe of Demand(Q) dQ - (Pe * Qe)
    integral_demand, _ = quad(demand_p, 0, Qe)
    consumer_surplus = integral_demand - (Pe * Qe)

    # 4. Calculate Producer Surplus (PS)
    # PS = (Pe * Qe) - Integral from q_min_supply to Qe of Supply(Q) dQ
    integral_supply, _ = quad(supply_p, q_min_supply, Qe)
    producer_surplus = (Pe * Qe) - integral_supply

    # 5. Calculate Total Welfare (TW)
    total_welfare = consumer_surplus + producer_surplus

    # Print the results in a step-by-step manner
    print("Finding the market equilibrium and total welfare:")
    print("-" * 50)
    print(f"Supply function: P = ln(Q^3 - 2)")
    print(f"Demand function: P = 18 * e^(-arctan(Q))")
    print("-" * 50)
    print("Step 1: Find Market Equilibrium (Qe, Pe)")
    print(f"Equilibrium Quantity (Qe) = {Qe:.4f}")
    print(f"Equilibrium Price (Pe) = {Pe:.4f}\n")
    
    print("Step 2: Calculate Consumer and Producer Surplus")
    print(f"Consumer Surplus (CS) = {consumer_surplus:.4f}")
    print(f"Producer Surplus (PS) = {producer_surplus:.4f}\n")

    print("Step 3: Calculate Total Welfare")
    # As requested, showing the numbers in the final equation
    print(f"Total Welfare = Consumer Surplus + Producer Surplus")
    print(f"Total Welfare = {consumer_surplus:.4f} + {producer_surplus:.4f} = {total_welfare:.4f}")

if __name__ == "__main__":
    main()

<<<24.7797>>>