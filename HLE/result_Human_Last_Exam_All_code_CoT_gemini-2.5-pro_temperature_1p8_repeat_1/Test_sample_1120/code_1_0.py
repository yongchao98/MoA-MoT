import numpy as np
from scipy.optimize import fsolve

def solve_and_print_excess_demand():
    """
    This function solves for the monopolist's optimal quantity and calculates the excess demand.
    """
    # Step 1: Define the monopolist's profit maximization problem.
    # The inverse demand curve is P(Q) = 4 + 0.03*Q^2 - 0.005*Q^3.
    # Total Revenue TR(Q) = P(Q) * Q = 4*Q + 0.03*Q^3 - 0.005*Q^4.
    # Marginal Cost MC = 0.
    # To maximize profit, we set Marginal Revenue (MR) equal to Marginal Cost (MC).
    # MR(Q) = d(TR)/dQ = 4 + 0.09*Q^2 - 0.02*Q^3.
    
    def marginal_revenue(q):
        """Calculates the marginal revenue for a given quantity q."""
        return 4 + 0.09 * q**2 - 0.02 * q**3

    # We need to solve MR(Q) = 0 to find the potential profit-maximizing quantity.
    # We use fsolve, a numerical solver, with an initial guess (e.g., 5).
    optimal_q_candidate_array = fsolve(marginal_revenue, 5)
    optimal_q_candidate = optimal_q_candidate_array[0]

    # Step 2: Verify this quantity is the profit-maximizing one within the [0, 10] constraint.
    def profit(q):
        """Calculates profit for a given quantity q."""
        if 0 <= q <= 10:
            return 4 * q + 0.03 * q**3 - 0.005 * q**4
        return -np.inf # Outside capacity, profit is effectively negative infinity

    profit_at_candidate = profit(optimal_q_candidate)
    profit_at_capacity = profit(10)

    # The monopolist chooses the quantity that gives higher profit.
    if profit_at_candidate > profit_at_capacity:
        quantity_demanded = optimal_q_candidate
    else:
        quantity_demanded = 10

    # Step 3: Calculate excess demand.
    # Quantity Demanded is the monopolist's chosen quantity.
    # Quantity Supplied is interpreted as the producer's maximum capacity.
    quantity_supplied = 10.0
    
    excess_demand = quantity_demanded - quantity_supplied

    # Output the final equation and the result
    print("The equilibrium is found where the monopolist maximizes profit. The quantity demanded at the monopolist's chosen price is the profit-maximizing quantity.")
    print("The quantity supplied is interpreted as the maximum capacity of the producer.")
    print("\nThe equation for excess demand is: Quantity Demanded - Quantity Supplied")
    print(f"{quantity_demanded:.4f} - {quantity_supplied:.4f} = {excess_demand:.4f}")

# Execute the function to find the answer
solve_and_print_excess_demand()