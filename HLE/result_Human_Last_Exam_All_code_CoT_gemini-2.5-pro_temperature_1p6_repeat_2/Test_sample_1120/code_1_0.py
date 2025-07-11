import numpy as np

def solve_market_equilibrium():
    """
    This function solves for the excess demand in the given economic scenario.
    """
    
    # Step 1 & 2: Define the relationship between Price (P) and Quantity Supplied (Q_s)
    # The demand for one customer is: q_i = 400 - 100*P + Q/100 + 3*Q^2 - Q^3/20
    # Let's assume the 'Q' that influences demand is the quantity supplied, Q_s.
    # Total demand Q_d = 100 * q_i = 40000 - 10000*P + Q_s + 300*Q_s^2 - 5*Q_s^3
    # For the producer to be able to sell Q_s, demand must be at least Q_s.
    # Q_d >= Q_s => 40000 - 10000*P + Q_s + 300*Q_s^2 - 5*Q_s^3 >= Q_s
    # This simplifies to 10000*P <= 40000 + 300*Q_s^2 - 5*Q_s^3
    # P <= 4 + 0.03*Q_s^2 - 0.0005*Q_s^3
    # To maximize profit, the producer chooses the highest possible price, making it an equality.
    # This gives us the inverse demand function P(Q_s).

    # Step 3: Set up and solve the producer's profit maximization problem.
    # The producer's profit is Pi = P * Q_s, since marginal cost is 0.
    # Pi(Q_s) = (4 + 0.03*Q_s^2 - 0.0005*Q_s^3) * Q_s
    # Pi(Q_s) = 4*Q_s + 0.03*Q_s^3 - 0.0005*Q_s^4
    # The producer can supply a maximum of 10 units, so Q_s is in [0, 10].
    # To find the maximum, we can check the derivative of Pi(Q_s) w.r.t Q_s (Marginal Revenue).
    # MR(Q_s) = 4 + 0.09*Q_s^2 - 0.002*Q_s^3
    # For Q_s in [0, 10], MR(Q_s) is always positive. For example, MR(0)=4 and MR(10)=11.
    # An increasing profit function on [0, 10] means the maximum is at the upper bound.
    
    Q_s_eq = 10.0
    print(f"Step 1: The producer maximizes profit by supplying the maximum possible quantity, which is Q_s = {Q_s_eq}")

    # Step 4: Calculate the equilibrium price based on the optimal quantity.
    P_eq = 4 + 0.03 * Q_s_eq**2 - 0.0005 * Q_s_eq**3
    print(f"Step 2: The equilibrium price set by the producer is P_eq = {P_eq}")

    # Step 5: Calculate the total quantity demanded at the equilibrium price and supplied quantity.
    Q_d_eq = 40000 - 10000 * P_eq + Q_s_eq + 300 * Q_s_eq**2 - 5 * Q_s_eq**3
    print(f"Step 3: At this price and supplied quantity, the total market demand is Q_d = {Q_d_eq}")

    # Step 6: Calculate and print the excess demand.
    excess_demand = Q_d_eq - Q_s_eq
    print(f"\nFinal Calculation:")
    print(f"Excess Demand = Quantity Demanded - Quantity Supplied")
    print(f"Excess Demand = {Q_d_eq} - {Q_s_eq}")
    print(f"The value of the excess demand is {excess_demand}.")

solve_market_equilibrium()