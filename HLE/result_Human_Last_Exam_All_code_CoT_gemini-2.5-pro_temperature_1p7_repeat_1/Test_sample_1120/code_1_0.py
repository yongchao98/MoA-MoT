import numpy as np

def solve_excess_demand():
    """
    Solves for the excess demand at the equilibrium price.
    """
    # Step 1: Define the inverse demand curve P(Q)
    # The individual demand is q_i = 400 - 100P + Q/100 + 3Q^2 - Q^3/20
    # Total demand Q = 100 * q_i.
    # Q = 100 * (400 - 100P + Q/100 + 3Q^2 - Q^3/20)
    # Q = 40000 - 10000P + Q + 300Q^2 - 5Q^3
    # 10000P = 40000 + 300Q^2 - 5Q^3
    # P(Q) = 4 + 0.03*Q^2 - 0.0005*Q^3
    def p_of_q(q):
        return 4 + 0.03 * q**2 - 0.0005 * q**3

    # Step 2 & 3: Find the producer's optimal quantity and the equilibrium price.
    # The producer maximizes Total Revenue TR(Q) = Q * P(Q) for Q in [0, 10].
    # TR(Q) = 4Q + 0.03*Q^3 - 0.0005*Q^4
    # Marginal Revenue MR(Q) = d(TR)/dQ = 4 + 0.09*Q^2 - 0.002*Q^3
    # On the interval [0, 10], MR(Q) is always positive, so TR(Q) is increasing.
    # Therefore, the producer chooses the maximum possible quantity.
    q_supplied = 10.0
    
    # The equilibrium price is the price to sell 10 units.
    p_equilibrium = p_of_q(q_supplied)

    print(f"The producer's optimal quantity supplied is Q_s = {q_supplied}")
    print(f"The equilibrium price set by the producer is P_e = {p_equilibrium}")

    # Step 4: Find all possible quantities demanded at the equilibrium price.
    # We need to solve P(Q) = p_equilibrium for Q.
    # 6.5 = 4 + 0.03*Q^2 - 0.0005*Q^3
    # This rearranges to the cubic equation: 0.0005*Q^3 - 0.03*Q^2 + 2.5 = 0
    # The coefficients for the polynomial Q^3, Q^2, Q^1, Q^0 are:
    coeffs = [0.0005, -0.03, 0, 2.5]
    roots = np.roots(coeffs)
    
    positive_real_roots = [r.real for r in roots if r.imag == 0 and r.real > 0]
    
    print(f"\nAt P = {p_equilibrium}, the possible levels of quantity demanded are: {', '.join(f'{r:.4f}' for r in sorted(positive_real_roots))}")

    # Step 5: Select the stable demand equilibrium.
    # Stability requires the demand curve to be downward sloping, i.e., dP/dQ < 0.
    # dP/dQ = 0.06*Q - 0.0015*Q^2
    def dp_dq(q):
        return 0.06 * q - 0.0015 * q**2

    stable_q_demanded = None
    for q in positive_real_roots:
        slope = dp_dq(q)
        print(f"At Q_d = {q:.4f}, the slope dP/dQ is {slope:.4f}.")
        if slope < 0:
            stable_q_demanded = q
            print(f"This corresponds to a stable demand equilibrium.")
        else:
            print(f"This corresponds to an unstable demand equilibrium.")

    if stable_q_demanded is None:
        print("No stable demand equilibrium found.")
        return
        
    print(f"\nThe stable quantity demanded at the equilibrium price is Q_d = {stable_q_demanded:.4f}")

    # Step 6: Calculate and print the excess demand.
    excess_demand = stable_q_demanded - q_supplied
    
    print("\nThe excess demand is the difference between the stable quantity demanded and the quantity supplied.")
    print("Final Equation:")
    print(f"{stable_q_demanded:.4f} - {q_supplied:.4f} = {excess_demand:.4f}")

    return excess_demand

# Execute the function to find the answer.
final_answer = solve_excess_demand()
print(f"\nThe final calculated value for excess demand is {final_answer:.4f}")
<<<48.5410>>>