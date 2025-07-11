import numpy as np
import math

def solve_excess_demand():
    """
    This function calculates the excess demand at the equilibrium price based on the problem's conditions.
    """
    # Step 1 & 2: Determine the producer's optimal quantity supplied (Q_s).
    # The producer's profit is π(Q) = P(Q) * Q. The inverse demand is P(Q) = 4 + 0.03*Q^2 - 0.0005*Q^3.
    # The profit function π(Q) = 4Q + 0.03Q^3 - 0.0005Q^4 is maximized on the interval [0, 10]
    # by choosing the largest possible quantity, because its derivative (marginal revenue)
    # is positive throughout this interval.
    Q_s = 10.0

    # Step 3: Calculate the equilibrium price (P_e) from the inverse demand function at Q_s.
    P_e = 4 + 0.03 * (Q_s**2) - 0.0005 * (Q_s**3)

    # Step 4: Determine the quantity demanded (Q_d) at P_e.
    # This requires solving the cubic demand equation for Q given P = P_e:
    # 100*P_e = 400 + 3*Q_d^2 - (Q_d^3)/20
    # Rearranging gives: Q_d^3 - 60*Q_d^2 + 5000 = 0
    # The coefficients for the cubic equation (ax^3 + bx^2 + cx + d = 0) are:
    coeffs = [1, -60, 0, 5000]
    roots = np.roots(coeffs)

    # Step 5: Identify the stable quantity demanded.
    # A stable equilibrium occurs on a downward-sloping demand curve.
    # For this specific demand function, this condition holds when Q > 40.
    # We find the root that is real, positive, and greater than 40.
    Q_d_stable = 0.0
    for r in roots:
        if np.isreal(r) and r > 40:
            Q_d_stable = float(np.real(r))
            break
            
    # Step 6: Calculate the excess demand.
    excess_demand = Q_d_stable - Q_s

    print("The final calculation for excess demand is based on the equation:")
    print("Excess Demand = Quantity Demanded (Q_d) - Quantity Supplied (Q_s)")
    print("-" * 50)
    print(f"The quantity supplied by the profit-maximizing producer is Q_s = {Q_s}")
    print(f"The stable quantity demanded at the equilibrium price is Q_d = {Q_d_stable:.4f}")
    print("\nThe final equation with values is:")
    # The request is to output each number in the final equation.
    print(f"{Q_d_stable:.4f} - {Q_s} = {excess_demand:.4f}")

solve_excess_demand()