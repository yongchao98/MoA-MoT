import numpy as np

def solve_excess_demand():
    """
    This function solves for the excess demand in the given economic scenario.
    """
    # Step 1: Define problem parameters and functions
    # The seller's maximum supply capacity
    quantity_supplied = 10

    # The inverse demand curve is derived from the individual demand:
    # q_i = 400 - 100P + Q/100 + 3Q^2 - Q^3/20
    # Q = 100 * q_i  => Q = 100 * (400 - 100P + Q/100 + 3Q^2 - Q^3/20)
    # Q = 40000 - 10000P + Q + 300Q^2 - 5Q^3
    # 10000P = 40000 + 300Q^2 - 5Q^3
    # P(Q) = 4 + 0.03*Q^2 - 0.0005*Q^3

    # Step 2: Determine the seller's optimal production quantity
    # Profit (pi) = Total Revenue (TR) since Total Cost is 0.
    # TR(Q) = P(Q) * Q = 4Q + 0.03Q^3 - 0.0005Q^4
    # Marginal Revenue MR(Q) = d(TR)/dQ = 4 + 0.09Q^2 - 0.002Q^3
    # To check if profit is increasing in the feasible range [0, 10], we evaluate MR.
    # The derivative of MR, d(MR)/dQ = 0.18Q - 0.006Q^2, is positive for Q in (0, 10].
    # Since MR(0) = 4 and MR(Q) is increasing on [0, 10], MR is always positive.
    # This means profit is maximized at the production capacity limit.
    print(f"The seller's profit is maximized by producing at full capacity.")
    print(f"Quantity Supplied (Q_s) = {quantity_supplied}")
    print("-" * 30)

    # Step 3: Determine the equilibrium price set by the seller
    # The seller sets the price to sell exactly 10 units.
    equilibrium_price = 4 + 0.03 * (quantity_supplied**2) - 0.0005 * (quantity_supplied**3)
    print(f"To sell {quantity_supplied} units, the seller sets the equilibrium price (P_eq) to: ${equilibrium_price:.2f}")
    print("-" * 30)

    # Step 4: Find the quantity demanded at the equilibrium price
    # We need to solve P_eq = 4 + 0.03*Q^2 - 0.0005*Q^3 for Q.
    # 6.5 = 4 + 0.03*Q^2 - 0.0005*Q^3
    # This simplifies to the cubic equation: Q^3 - 60Q^2 + 5000 = 0
    # The coefficients for the polynomial Q^3 - 60Q^2 + 0Q + 5000 are [1, -60, 0, 5000].
    coeffs = [1, -60, 0, 5000]
    roots = np.roots(coeffs)
    print(f"At P_eq = ${equilibrium_price:.2f}, we solve for the quantity demanded (Q_d) from the equation:")
    print("Q_d^3 - 60*Q_d^2 + 5000 = 0")
    
    positive_real_roots = [np.real(r) for r in roots if np.isreal(r) and np.real(r) > 0]
    print(f"The positive real roots are approximately: {', '.join([f'{r:.2f}' for r in sorted(positive_real_roots)])}")
    print("-" * 30)

    # Step 5: Identify the stable market equilibrium
    # The stable equilibrium is where the demand curve is downward sloping (dP/dQ < 0).
    # dP/dQ = 0.06*Q - 0.0015*Q^2
    quantity_demanded = 0
    print("Checking stability for each positive root (stable if dP/dQ < 0):")
    for q in sorted(positive_real_roots):
        slope = 0.06 * q - 0.0015 * (q**2)
        stability = "stable" if slope < 0 else "unstable"
        print(f"For Q = {q:.2f}, dP/dQ = {slope:.2f} -> ({stability})")
        if slope < 0:
            quantity_demanded = q
    
    print(f"\nThe market will coordinate on the stable equilibrium.")
    print(f"Quantity Demanded (Q_d) = {quantity_demanded:.2f}")
    print("-" * 30)

    # Step 6: Calculate and print the final excess demand
    excess_demand = quantity_demanded - quantity_supplied
    print("Calculating the excess demand:")
    print(f"Excess Demand = Quantity Demanded - Quantity Supplied")
    print(f"Excess Demand = {quantity_demanded:.2f} - {quantity_supplied} = {excess_demand:.2f}")

solve_excess_demand()