import numpy as np

def solve_excess_demand():
    """
    This function calculates the excess demand by following these steps:
    1. Determines the profit-maximizing quantity for the producer to supply.
    2. Calculates the price the producer will set to sell that quantity (the equilibrium price).
    3. Solves for the quantity consumers will demand at that price, considering equilibrium stability.
    4. Computes the excess demand as the difference between quantity demanded and supplied.
    """

    # Step 1 & 2: Determine the producer's optimal quantity and the resulting equilibrium price.
    # The producer's profit is pi(Q) = P(Q) * Q. With zero marginal cost and a capacity of 10,
    # the producer's goal is to maximize Total Revenue, TR(Q), for Q <= 10.
    # The inverse demand curve is P(Q) = 4 + 0.03*Q**2 - 0.0005*Q**3.
    # Total Revenue is TR(Q) = 4*Q + 0.03*Q**3 - 0.0005*Q**4.
    # Marginal Revenue is MR(Q) = 4 + 0.09*Q**2 - 0.002*Q**3.
    # For any quantity Q in [0, 10], MR(Q) is positive (e.g., MR(0)=4, MR(10)=11).
    # Since MR > MC (which is 0), profit increases with quantity.
    # Therefore, the producer supplies its maximum capacity.
    q_supply = 10.0

    # The producer sets the price to sell exactly 10 units.
    # This is the equilibrium price P_eq = P(10).
    p_eq = 4 + 0.03 * q_supply**2 - 0.0005 * q_supply**3

    print(f"The producer will supply Q_supply = {q_supply:.1f} units.")
    print(f"To sell this quantity, the producer sets an equilibrium price P_eq = {p_eq:.2f}.")

    # Step 3: Find the quantity demanded by the market at P_eq.
    # We need to solve for Q in the equation P(Q) = p_eq.
    # 6.5 = 4 + 0.03*Q**2 - 0.0005*Q**3
    # This rearranges to the cubic equation: Q**3 - 60*Q**2 + 5000 = 0.
    coeffs = [1, -60, 0, 5000]
    roots = np.roots(coeffs)

    print(f"\nAt the price of {p_eq:.2f}, the possible quantities demanded are the roots of the equation Q³ - 60Q² + 5000 = 0.")
    print(f"The roots are: {[f'{r:.2f}' for r in sorted(roots)]}.")

    # Due to the nature of the demand curve, we must choose the stable equilibrium.
    # For this model, an equilibrium is stable if Q > 40.
    # We filter out negative and unstable roots.
    stable_roots = [r.real for r in roots if np.isreal(r) and r.real > 40]
    q_demand = stable_roots[0]

    print(f"The stable quantity demanded, Q_demand, is {q_demand:.2f}.")

    # Step 4: Calculate the excess demand.
    excess_demand = q_demand - q_supply

    print("\nExcess demand is the difference between the stable quantity demanded and the quantity supplied.")
    # The final equation output requested by the user
    print(f"Final Calculation: {q_demand:.2f} - {q_supply:.1f} = {excess_demand:.2f}")

    return excess_demand

if __name__ == '__main__':
    solve_excess_demand()