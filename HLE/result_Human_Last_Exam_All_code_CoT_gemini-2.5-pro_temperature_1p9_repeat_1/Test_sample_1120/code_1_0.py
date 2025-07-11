import numpy as np
import math

def solve_excess_demand():
    """
    This function calculates the excess demand by following these steps:
    1. Determines the producer's optimal supply quantity.
    2. Calculates the price the producer will set.
    3. Finds the stable market demand at that price.
    4. Computes the excess demand.
    """

    # Step 1 & 2: Determine the producer's profit-maximizing supply quantity.
    # The relationship between P and Q is derived from the demand function:
    # q_i = 400 - 100*P + Q/100 + 3*Q^2 - Q^3/20
    # Since Q = 100*q_i, we substitute q_i = Q/100 into the equation.
    # This leads to the inverse demand curve: P(Q) = 4 + 0.03*Q^2 - 0.0005*Q^3.
    #
    # The producer's profit is Total Revenue (TR) minus Total Cost (TC).
    # Since marginal cost is 0 for Q <= 10, TC=0 in the feasible range.
    # So, the producer maximizes TR = P(Q) * Q = 4Q + 0.03Q^3 - 0.0005Q^4.
    #
    # The marginal revenue is MR(Q) = d(TR)/dQ = 4 + 0.09Q^2 - 0.002Q^3.
    # By analyzing MR(Q) on the interval [0, 10], we find it is always positive.
    # This means TR(Q) is always increasing on [0, 10].
    # Therefore, the producer maximizes revenue by supplying the maximum possible quantity.
    Q_supply = 10.0
    print(f"Step 1: The producer's optimal quantity to supply is Q_supply = {Q_supply}")

    # Step 3: Calculate the equilibrium price set by the producer.
    # The price is set using the inverse demand curve at Q_supply = 10.
    P_eq = 4 + 0.03 * (Q_supply**2) - 0.0005 * (Q_supply**3)
    print(f"Step 2: The producer sets the equilibrium price P_eq = {P_eq}")

    # Step 4: Determine the actual quantity demanded at the equilibrium price.
    # At P_eq = 6.5, we solve for Q in the demand relationship:
    # 0 = 400 - 100*P_eq + 3*Q^2 - Q^3/20
    # 0 = 400 - 650 + 3*Q^2 - 0.05*Q^3
    # Q^3 - 60*Q^2 + 5000 = 0
    # We find the roots of this polynomial.
    coefficients = [1, -60, 0, 5000]
    roots = np.roots(coefficients)

    print(f"\nStep 3: Finding the market demand at P = {P_eq}.")
    print(f"The potential demand quantities (roots of the demand equation) are: {roots}")

    # Due to the externality, we must check which equilibrium is stable.
    # An equilibrium Q* is stable if G'(Q*) < 0, where G'(Q) = 600*Q - 15*Q^2.
    Q_demand = 0
    print("Checking stability of the positive real roots:")
    for r in roots:
        # We only consider real, positive roots for economic sense.
        if r.imag == 0 and r.real > 0:
            q_val = r.real
            # Stability check: G'(Q) = 15Q(40-Q). Stable if G'(Q) < 0.
            stability_metric = 15 * q_val * (40 - q_val)
            if stability_metric < 0:
                stability = 'stable'
                Q_demand = q_val
            else:
                stability = 'unstable'
            print(f"  - For Q = {q_val:.4f}, the equilibrium is {stability}.")
    
    if Q_demand == 0:
         print("No stable positive demand equilibrium found.")
         return

    print(f"The stable quantity demanded is Q_demand = {Q_demand}")

    # Step 5: Calculate and print the excess demand.
    excess_demand = Q_demand - Q_supply
    print("\nStep 4: The final equation for excess demand is:")
    print(f"Excess Demand = Q_demand - Q_supply")
    print(f"Excess Demand = {Q_demand} - {Q_supply}")
    print(f"The value of the excess demand is: {excess_demand}")

    # The precise analytical value for the root is 25 + 15*sqrt(5)
    # The precise excess demand is 15 + 15*sqrt(5)
    precise_excess_demand = 15 + 15 * math.sqrt(5)
    return precise_excess_demand

if __name__ == '__main__':
    result = solve_excess_demand()
