import numpy as np

def solve_excess_demand():
    """
    This function solves the problem step-by-step as outlined in the plan.
    """
    # Step 1: Define the demand relationship derived from the problem description.
    # The individual demand is q_i = 400 - 100P + Q/100 + 3Q^2 - Q^3/20.
    # For 100 customers, total demand Q = 100*q_i, which leads to the
    # market relationship (inverse demand curve) P(Q):
    # 10000P = 40000 + 300Q^2 - 5Q^3  or P(Q) = 4 + 0.03*Q^2 - 0.0005*Q^3

    # Step 2: Determine the producer's optimal quantity.
    # Profit Pi(Q) = P(Q)*Q = 4Q + 0.03*Q^3 - 0.0005*Q^4 (since Cost=0 for Q<=10).
    # Marginal Revenue MR(Q) = d(Pi)/dQ = 4 + 0.09*Q^2 - 0.002*Q^3.
    # Marginal Cost MC = 0 for Q <= 10.
    # Let's check MR at the production limit Q=10.
    # MR(10) = 4 + 0.09*(10^2) - 0.002*(10^3) = 4 + 9 - 2 = 11.
    # Since MR(Q) is positive for all Q in [0, 10] and greater than MC=0,
    # the producer maximizes profit by producing the maximum possible quantity.
    q_supplied = 10.0
    print(f"The producer's optimal quantity to supply is {q_supplied} units.")

    # Step 3: Find the equilibrium price set by the producer.
    # P_eq = P(q_supplied)
    p_eq = 4 + 0.03 * q_supplied**2 - 0.0005 * q_supplied**3
    print(f"To sell this quantity, the producer sets the equilibrium price at P = ${p_eq:.2f}.")

    # Step 4 & 5: Find all possible demands at this price and select the stable one.
    # At P_eq, find all Q_d such that P(Q_d) = p_eq.
    # p_eq = 4 + 0.03*Q_d^2 - 0.0005*Q_d^3
    # This rearranges to the cubic equation: Q_d^3 - 60*Q_d^2 + 5000 = 0.
    coeffs = [1, -60, 0, 5000]
    roots = np.roots(coeffs)
    
    # The roots are the possible quantities demanded. We consider positive real roots.
    # The demand equilibrium is stable where the P(Q) curve is downward-sloping,
    # which means dP/dQ < 0.
    # dP/dQ = 0.06Q - 0.0015Q^2. This is negative for Q > 40.
    # One root is 10 (unstable, as dP/dQ(10) > 0). The largest root will be the stable one.
    q_demanded_stable = max(r.real for r in roots if r.imag == 0 and r.real > 0)
    print(f"At this price, due to network effects, the stable market demand is {q_demanded_stable:.4f} units.")
    
    # Step 6: Calculate and print the excess demand.
    excess_demand = q_demanded_stable - q_supplied
    
    print("\nExcess demand is the stable quantity demanded minus the quantity supplied.")
    print("The final equation for excess demand is:")
    print(f"{q_demanded_stable:.4f} - {q_supplied:.4f} = {excess_demand:.4f}")

solve_excess_demand()