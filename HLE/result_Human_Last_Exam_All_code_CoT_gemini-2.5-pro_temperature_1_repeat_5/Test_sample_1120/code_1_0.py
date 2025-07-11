import numpy as np

def solve_excess_demand():
    """
    Calculates the excess demand at the equilibrium price based on the problem's conditions.
    """
    # Step 1: Define producer's profit-maximizing quantity (Q_s)
    # The producer's profit is TR(Q) = P(Q)*Q = (4 + 0.03Q^2 - 0.0005Q^3)*Q.
    # The marginal revenue is MR(Q) = 4 + 0.09Q^2 - 0.002Q^3.
    # For Q in [0, 10], MR(Q) is always positive, so revenue is maximized at the capacity limit.
    Q_s = 10.0

    # Step 2: Find the equilibrium price (P_eq) set by the producer to sell Q_s
    # P(Q) = 4 + 0.03*Q^2 - 0.0005*Q^3
    P_eq = 4 + 0.03 * (Q_s**2) - 0.0005 * (Q_s**3)

    # Step 3: Find all possible quantities demanded (Q) at P_eq
    # The aggregate demand relation is 10000*P = 40000 + 300*Q^2 - 5*Q^3.
    # Substituting P_eq = 6.5 gives: 65000 = 40000 + 300*Q^2 - 5*Q^3
    # Rearranging gives a cubic equation: Q^3 - 60*Q^2 + 5000 = 0.
    # Coefficients for the cubic equation Q^3 - 60*Q^2 + 0*Q + 5000 = 0
    coeffs = [1, -60, 0, 5000]
    roots = np.roots(coeffs)

    # Filter for real, positive roots which are economically meaningful
    possible_Qd = [r.real for r in roots if np.isreal(r) and r.real >= 0]

    # Step 4: Identify the stable demand equilibrium (Q_d)
    # The equilibrium is stable if the slope of the inverse demand curve, dP/dQ, is negative.
    # dP/dQ = 0.06*Q - 0.0015*Q^2
    stable_Qd = 0
    for q in possible_Qd:
        # Check stability condition
        slope = 0.06 * q - 0.0015 * (q**2)
        if slope < 0:
            stable_Qd = q
            break # Assume only one stable equilibrium exists

    # Step 5: Calculate and print the excess demand
    if stable_Qd > 0:
        excess_demand = stable_Qd - Q_s
        print(f"The producer supplies Q_s = {Q_s:.2f} units.")
        print(f"The equilibrium price set by the producer is P_eq = {P_eq:.2f}.")
        print(f"At this price, the stable quantity demanded by the market is Q_d = {stable_Qd:.2f}.")
        print("\nThe excess demand is the difference between the stable quantity demanded and the quantity supplied.")
        print(f"Excess Demand = Q_d - Q_s")
        print(f"Final Equation: {stable_Qd:.2f} - {Q_s:.2f} = {excess_demand:.2f}")
    else:
        print("No stable demand equilibrium found.")

solve_excess_demand()