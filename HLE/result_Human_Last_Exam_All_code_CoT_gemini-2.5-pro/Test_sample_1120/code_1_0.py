import numpy as np

def solve_excess_demand():
    """
    This function calculates the excess demand by following these steps:
    1. Determines the profit-maximizing quantity for the producer to supply.
    2. Calculates the equilibrium price based on that quantity.
    3. Finds the total quantity demanded by consumers at that price.
    4. Computes the excess demand as the difference between quantity demanded and supplied.
    """

    # Step 1: Determine the Producer's Supply (Q_s)
    # The inverse demand curve is P(Q) = 4 + 0.03*Q^2 - 0.0005*Q^3.
    # Total Revenue TR(Q) = P(Q) * Q = 4*Q + 0.03*Q^3 - 0.0005*Q^4.
    # Marginal Revenue MR(Q) = d(TR)/dQ = 4 + 0.09*Q^2 - 0.002*Q^3.
    # The producer's marginal cost is 0 for Q <= 10.
    # To maximize profit, the producer chooses Q where MR(Q) = MC = 0, or at the
    # capacity limit if MR > 0.
    # Let's check MR(Q) in the range [0, 10].
    # MR(10) = 4 + 0.09*(10**2) - 0.002*(10**3) = 4 + 9 - 2 = 11.
    # Since MR(Q) is positive for all Q in [0, 10], the producer's revenue is
    # always increasing in this range. Therefore, to maximize profit, the
    # producer will supply the maximum possible quantity.
    Q_supplied = 10.0

    # Step 2: Find the Equilibrium Price (P_eq)
    # The producer sets the price the market will bear for Q_supplied = 10.
    P_eq = 4 + 0.03 * (Q_supplied**2) - 0.0005 * (Q_supplied**3)

    # Step 3: Calculate the Quantity Demanded (Q_d) at P_eq
    # We need to solve P(Q) = P_eq for Q.
    # P_eq = 4 + 0.03*Q^2 - 0.0005*Q^3
    # This can be rearranged into the cubic equation:
    # 0.0005*Q^3 - 0.03*Q^2 + (P_eq - 4) = 0
    # Or, multiplying by 2000: Q^3 - 60*Q^2 + 5000 = 0.
    coeffs = [1, -60, 0, 5000]
    roots = np.roots(coeffs)

    # We are interested in real, non-negative solutions for quantity.
    # One root is Q=10 (the quantity supplied), but there might be others.
    # The derivative of the inverse demand curve, dP/dQ = 0.06Q - 0.0015Q^2,
    # indicates that the demand curve is not always downward sloping. Stable
    # market equilibria occur on the downward-sloping portions.
    # The largest positive real root corresponds to the stable equilibrium.
    Q_demanded = max(r.real for r in roots if np.isreal(r) and r.real >= 0)

    # Step 4: Calculate and Print the Excess Demand
    excess_demand = Q_demanded - Q_supplied

    print(f"The producer will supply Q_supplied = {Q_supplied:.2f} units.")
    print(f"The resulting equilibrium price is P_eq = {P_eq:.2f}.")
    print(f"At this price, the stable quantity demanded is Q_demanded = {Q_demanded:.2f}.")
    print("\nThe excess demand is calculated as:")
    print(f"Excess Demand = Quantity Demanded - Quantity Supplied")
    print(f"Excess Demand = {Q_demanded:.2f} - {Q_supplied:.2f}")
    print(f"\nThe final value of the excess demand is: {excess_demand:.2f}")

solve_excess_demand()
<<<48.54>>>