import numpy as np

def solve_excess_demand():
    """
    This function solves for the excess demand based on the problem description.
    """
    # Step 1 & 2: Define profit function and find optimal quantity.
    # The inverse demand curve is P(Q) = 4 + 0.03*Q**2 - 0.0005*Q**3.
    # The profit function is pi(Q) = P(Q) * Q = 4*Q + 0.03*Q**3 - 0.0005*Q**4.
    # The marginal revenue is MR(Q) = 4 + 0.09*Q**2 - 0.002*Q**3.
    # The marginal cost is MC = 0 for Q <= 10.
    # The producer maximizes profit by producing as long as MR > MC.
    # Let's check MR at the capacity limit Q=10.
    Q_check = 10
    MR_at_10 = 4 + 0.09 * Q_check**2 - 0.002 * Q_check**3
    
    # Since MR(10) = 11, which is greater than MC=0, and MR is increasing on [0, 10],
    # the producer's optimal strategy is to produce at full capacity.
    Q_supply = 10
    print(f"The producer's profit-maximizing quantity to supply is Q_supply = {Q_supply}")

    # Step 3: Determine the equilibrium price for Q_supply = 10.
    # P = 4 + 0.03*Q**2 - 0.0005*Q**3
    P_eq = 4 + 0.03 * Q_supply**2 - 0.0005 * Q_supply**3
    print(f"The equilibrium price set by the producer is P = {P_eq}")

    # Step 4: Determine the quantity demanded at the equilibrium price.
    # The quantity demanded Q_d must satisfy: 10000*P = 40000 + 300*Q**2 - 5*Q**3
    # 10000 * 6.5 = 40000 + 300*Q_d**2 - 5*Q_d**3
    # This gives the cubic equation: Q_d**3 - 60*Q_d**2 + 5000 = 0
    # The coefficients for the polynomial are [1, -60, 0, 5000].
    coeffs = [1, -60, 0, 5000]
    roots = np.roots(coeffs)
    print(f"At price P={P_eq}, the possible consistent market demand quantities are: {np.round(roots, 2)}")

    # The producer's action of supplying 10 units creates a focal point for consumers'
    # expectations. The rational market outcome is the one where demand matches this supply.
    Q_demand = 10.0
    print(f"The relevant quantity demanded, given the producer's action, is Q_demand = {Q_demand}")

    # Step 5: Calculate the excess demand.
    # Excess Demand = Quantity Demanded - Quantity Supplied
    excess_demand = Q_demand - Q_supply
    
    print("\n--- Calculation of Excess Demand ---")
    print(f"The final equation for excess demand is: {Q_demand} - {Q_supply} = {excess_demand}")
    print(f"The value of the excess demand at the equilibrium price is {excess_demand}")

solve_excess_demand()