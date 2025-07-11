import sympy

def solve_excess_demand():
    """
    This function solves for the excess demand at the equilibrium price by following these steps:
    1. Derives the market inverse demand curve P(Q).
    2. Defines the producer's profit function Profit(Q) = P(Q) * Q.
    3. Finds the profit-maximizing quantity for the producer, considering the Q<=10 constraint.
    4. Determines the equilibrium quantity and price.
    5. Calculates the excess demand, which is Quantity Demanded - Quantity Supplied.
    """
    # Define symbols for price and quantity
    P, Q = sympy.symbols('P Q')

    # Step 1: Derive the market inverse demand curve.
    # The individual demand is q_i = 400 - 100*P + Q/100 + 3*Q**2 - Q**3/20.
    # With 100 identical customers, total demand Q = 100 * q_i, which means q_i = Q/100.
    # We set up the equation by substituting q_i:
    # Q/100 = 400 - 100*P + Q/100 + 3*Q**2 - Q**3/20
    # The Q/100 terms cancel, and we rearrange to solve for P in terms of Q.
    # 100*P = 400 + 3*Q**2 - Q**3/20
    P_demand = (400 + 3*Q**2 - Q**3/20) / 100
    
    print("Step 1: The derived market inverse demand curve P(Q) is:")
    sympy.pprint(P_demand)
    print("-" * 40)

    # Step 2: Define the producer's profit function.
    # With a marginal cost of 0, Profit = Total Revenue = P * Q.
    profit = P_demand * Q
    
    print("Step 2: The producer's profit function Profit(Q) is:")
    sympy.pprint(sympy.expand(profit))
    print("-" * 40)

    # Step 3: Find the profit-maximizing quantity.
    # We take the derivative of the profit function with respect to Q.
    profit_derivative = sympy.diff(profit, Q)
    
    print("Step 3: The derivative of the profit function, d(Profit)/dQ, is:")
    sympy.pprint(profit_derivative)
    
    # The producer can only sell up to Q=10. We check the slope of the profit function at this upper limit.
    derivative_at_10 = profit_derivative.subs(Q, 10)
    print(f"\nAt the maximum capacity Q=10, the value of the derivative is: {derivative_at_10}")
    print("Since this value is positive, the producer's profit is still increasing at Q=10.")
    print("Therefore, to maximize profit, the producer will sell the maximum possible quantity.")
    print("-" * 40)

    # Step 4: Determine the equilibrium quantity and price.
    # The quantity supplied by the producer is their maximum capacity.
    Q_s = 10
    print(f"Step 4: The equilibrium quantity supplied by the producer is Q_s = {Q_s}.")

    # The equilibrium price is the price that makes consumers demand exactly Q_s=10 units.
    P_eq = P_demand.subs(Q, Q_s)
    print(f"The equilibrium price is P_eq = P({Q_s}) = {float(P_eq)}.")
    print("-" * 40)

    # Step 5: Calculate the excess demand.
    # At the equilibrium price P_eq, the quantity demanded Q_d is, by definition of the demand curve, equal to Q_s.
    Q_d = Q_s
    print(f"Step 5: At P_eq = {float(P_eq)}, the quantity demanded is Q_d = {Q_d}.")
    
    excess_demand = Q_d - Q_s
    print("Excess demand is the difference between quantity demanded and quantity supplied.")
    print("\nFinal Calculation:")
    print(f"{Q_d} - {Q_s} = {excess_demand}")

solve_excess_demand()