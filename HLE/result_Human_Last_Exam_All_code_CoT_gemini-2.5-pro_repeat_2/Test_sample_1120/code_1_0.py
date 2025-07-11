import numpy as np

def solve_excess_demand():
    """
    Calculates the excess demand at the equilibrium price based on the problem description.
    """
    # Step 1: Derive the Inverse Market Demand Curve.
    # The demand for one customer is q_i = 400 - 100P + Q/100 + 3Q^2 - Q^3/20.
    # Total market demand Q = 100 * q_i.
    # Substituting q_i = Q/100:
    # Q/100 = 400 - 100P + Q/100 + 3Q^2 - Q^3/20
    # 0 = 400 - 100P + 3Q^2 - Q^3/20
    # The above equation is for a single customer based on total market Q.
    # Let's use the total demand equation directly:
    # Q_d = 100 * (400 - 100P + Q_d/100 + 3*(Q_d)^2 - (Q_d)^3/20)
    # Q_d = 40000 - 10000P + Q_d + 300*(Q_d)^2 - 5*(Q_d)^3
    # 0 = 40000 - 10000P + 300*(Q_d)^2 - 5*(Q_d)^3
    # 10000P = 40000 + 300*(Q_d)^2 - 5*(Q_d)^3
    # P(Q_d) = 4 + 0.03*(Q_d)^2 - 0.0005*(Q_d)^3
    print("Step 1: The inverse market demand curve is P(Q) = 4 + 0.03*Q^2 - 0.0005*Q^3.")

    # Step 2: Determine the Producer's Optimal Supply Quantity.
    # Profit (pi) = Revenue - Cost. Since Cost=0 for Q<=10, pi = Revenue = P(Q) * Q.
    # pi(Q) = (4 + 0.03*Q^2 - 0.0005*Q^3) * Q = 4Q + 0.03*Q^3 - 0.0005*Q^4.
    # To maximize profit, we check the marginal revenue MR = d(pi)/dQ.
    # MR(Q) = 4 + 0.09*Q^2 - 0.002*Q^3.
    # The marginal cost MC is 0 for Q <= 10.
    # We need to find where profit is maximized in the interval [0, 10].
    # Let's evaluate MR at the endpoint Q=10:
    # MR(10) = 4 + 0.09*(10^2) - 0.002*(10^3) = 4 + 9 - 2 = 11.
    # Since MR(Q) is positive for all Q in [0, 10], the profit function is always increasing on this interval.
    # Therefore, the producer maximizes profit by producing the maximum possible quantity.
    quantity_supplied = 10
    print(f"Step 2: The producer's optimal quantity to supply is Q_s = {quantity_supplied} units.")

    # Step 3: Find the Equilibrium Price and Quantity Demanded.
    # The equilibrium price is the price that makes the market demand equal to the quantity supplied.
    # We use the inverse demand curve with Q = 10.
    price_eq = 4 + 0.03 * (quantity_supplied**2) - 0.0005 * (quantity_supplied**3)
    print(f"Step 3: The equilibrium price P_eq is calculated at Q = {quantity_supplied}.")
    print(f"P_eq = 4 + 0.03*({quantity_supplied}^2) - 0.0005*({quantity_supplied}^3) = {price_eq}")
    
    # At this price, by definition, the quantity demanded by the market is equal to the quantity supplied.
    quantity_demanded = quantity_supplied
    print(f"At this price, the quantity demanded by the market is Q_d = {quantity_demanded} units.")

    # Step 4: Calculate Excess Demand.
    # Excess Demand = Quantity Demanded - Quantity Supplied
    excess_demand = quantity_demanded - quantity_supplied
    print("\nStep 4: Calculate the final excess demand.")
    print(f"Excess Demand = Quantity Demanded ({quantity_demanded}) - Quantity Supplied ({quantity_supplied}) = {excess_demand}")

solve_excess_demand()