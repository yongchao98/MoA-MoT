import numpy as np

def solve_excess_demand():
    """
    Calculates the excess demand at the equilibrium price for the given market conditions.
    """
    # Step 1: Model the producer's problem.
    # The producer's capacity is 10 units. We need to find the profit-maximizing quantity Q in [0, 10].
    # Profit = P(Q) * Q.

    # Step 2: Determine the market inverse demand curve P(Q).
    # Individual demand: q_i = 400 - 100P + Q/100 + 3Q^2 - Q^3/20
    # Total demand Q_d = 100 * q_i.
    # In equilibrium, the total quantity Q in the formula is the quantity demanded, so Q_d = Q.
    # Q = 100 * (400 - 100P + Q/100 + 3*Q**2 - Q**3/20)
    # Q = 40000 - 10000P + Q + 300*Q**2 - 5*Q**3
    # 0 = 40000 - 10000P + 300*Q**2 - 5*Q**3
    # 10000P = 40000 + 300*Q**2 - 5*Q**3
    # P(Q) = 4 + 0.03*Q**2 - 0.0005*Q**3

    def inverse_demand_p(Q):
        return 4 + 0.03 * Q**2 - 0.0005 * Q**3

    # Step 3: Find the profit-maximizing quantity Q.
    # Profit(Q) = P(Q) * Q = 4Q + 0.03*Q**3 - 0.0005*Q**4
    def profit(Q):
        return inverse_demand_p(Q) * Q

    # To maximize profit on the interval [0, 10], we check the endpoints.
    # The derivative of the profit function is the Marginal Revenue (MR).
    # MR(Q) = 4 + 0.09*Q^2 - 0.002*Q^3.
    # MR(Q) > 0 for all Q in [0, 10], meaning profit is always increasing.
    # Therefore, the maximum profit occurs at the capacity limit Q = 10.
    
    # Step 4: Find the equilibrium price and quantity.
    # The producer chooses to supply the profit-maximizing quantity.
    Q_s = 10
    print(f"The profit-maximizing quantity for the producer to supply is Q_s = {Q_s}.")

    # The equilibrium price is the price from the demand curve at Q = 10.
    P_e = inverse_demand_p(Q_s)
    print(f"To sell {Q_s} units, the producer must set the equilibrium price P_e = {P_e}.")

    # Step 5: Calculate the quantity demanded at the equilibrium price.
    # We use the original demand function to find the total quantity demanded (Q_d)
    # when P = P_e and the total market quantity Q is the equilibrium quantity Q_s.
    
    # Individual demand q_i
    q_i_demand = 400 - 100 * P_e + Q_s / 100 + 3 * Q_s**2 - Q_s**3 / 20
    
    # Total demand Q_d
    Q_d = 100 * q_i_demand
    print(f"At the equilibrium price of {P_e}, the total quantity demanded is Q_d = {Q_d}.")

    # Finally, calculate the excess demand.
    excess_demand = Q_d - Q_s

    print("\nThe final equation is: Excess Demand = Quantity Demanded - Quantity Supplied")
    print(f"Excess Demand = {Q_d} - {Q_s}")
    print(f"The value of the excess demand is: {excess_demand}")
    
    return excess_demand

if __name__ == '__main__':
    result = solve_excess_demand()
    # The final answer is enclosed in <<<>>>
    print(f"\n<<<{result}>>>")