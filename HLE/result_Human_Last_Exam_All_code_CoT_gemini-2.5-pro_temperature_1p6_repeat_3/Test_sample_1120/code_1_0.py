import numpy as np
from scipy.optimize import fsolve

def solve_excess_demand():
    """
    Calculates the excess demand at the equilibrium price for the given market conditions.
    """
    # Step 1 & 2: Define the inverse demand curve P(Q) and profit function.
    # From the problem statement: q_i = 400 - 100P + Q/100 + 3Q^2 - Q^3/20
    # With Q = 100 * q_i, we get Q/100 = 400 - 100P + Q/100 + 3Q^2 - Q^3/20
    # Simplifying gives: 100P = 400 + 3Q^2 - Q^3/20
    # So, the inverse demand curve P(Q) is:
    def inverse_demand(Q):
        return 4 + 0.03 * Q**2 - 0.005 * Q**3

    # Profit (pi) = Revenue - Cost. Since MC = 0, Profit = Revenue = P(Q) * Q
    def profit(Q):
        return inverse_demand(Q) * Q

    # Step 3: Find the profit-maximizing quantity Q_m by finding where d(pi)/dQ = 0
    # d(pi)/dQ = 4 + 0.09*Q^2 - 0.02*Q^3
    def profit_derivative(Q):
        return 4 + 0.09 * Q**2 - 0.02 * Q**3

    # Use a numerical solver to find the root of the profit derivative.
    # An initial guess of 7.5 is reasonable based on preliminary analysis.
    initial_guess = 7.5
    # fsolve returns an array, so we take the first element.
    optimal_quantity = fsolve(profit_derivative, initial_guess)[0]

    # The producer can only sell up to 10 units. We must check if the optimal
    # quantity is within this constraint and yields a higher profit than selling 10 units.
    # The derivative is positive before the root and negative after, indicating a maximum.
    # The profit at the capacity limit (Q=10) is profit(10) = 20.
    # The profit at the optimal quantity is profit(optimal_quantity).
    # Since profit(optimal_quantity) is higher, the producer chooses optimal_quantity.
    
    # In this equilibrium, the quantity the producer chooses to supply is this optimal quantity.
    quantity_supplied = optimal_quantity
    
    # Step 4: Determine the equilibrium price P_m.
    # This is the price that makes consumers demand exactly the quantity supplied.
    equilibrium_price = inverse_demand(quantity_supplied)

    # Step 5: Calculate the quantity demanded at the equilibrium price.
    # By definition, at the price the monopolist sets, the quantity demanded
    # is exactly the quantity they are supplying to maximize profit.
    quantity_demanded = quantity_supplied

    # Step 6: Calculate and print the excess demand.
    excess_demand = quantity_demanded - quantity_supplied

    print("--- Monopoly Equilibrium Calculation ---")
    print(f"Profit-Maximizing Quantity (Supply): Qs = {quantity_supplied:.4f}")
    print(f"Resulting Equilibrium Price: P = {equilibrium_price:.4f}")
    print(f"Quantity Demanded at this Price: Qd = {quantity_demanded:.4f}")
    print("\n--- Excess Demand Calculation ---")
    print(f"Excess Demand = Quantity Demanded - Quantity Supplied")
    # Final output as requested
    print(f"The final equation is: {quantity_demanded:.4f} - {quantity_supplied:.4f} = {excess_demand:.4f}")

solve_excess_demand()