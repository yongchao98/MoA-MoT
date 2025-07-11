import sympy

def demonstrate_contraction_proof():
    """
    This function demonstrates the proof that the Bellman operator is a contraction mapping,
    showing that the reward function does not affect the contraction property.
    """
    
    # Define symbolic variables
    # V1 and V2 are two different value functions
    V1, V2 = sympy.symbols('V1 V2', cls=sympy.Function)
    # T is the Bellman operator
    T = sympy.symbols('T', cls=sympy.Function)
    # s is a state
    s = sympy.symbols('s')
    # a is an action, a_star is the optimal action for TV1(s)
    a, a_star = sympy.symbols('a a_star')
    # R is the reward function
    R = sympy.symbols('R', cls=sympy.Function)
    # P is the transition probability function
    P_V1 = sympy.Symbol("sum[P(s'|s,a_star)V1(s')]")
    P_V2 = sympy.Symbol("sum[P(s'|s,a_star)V2(s')]")
    # gamma is the discount factor
    gamma = sympy.symbols('gamma')
    
    # The Bellman operator applied to V1 at state s
    # a_star is the action that maximizes the expression for V1
    TV1_s = R(s, a_star) + gamma * P_V1

    # The Bellman operator applied to V2 at state s
    # For the inequality, we use the same action a_star. The max for V2 might be at another action.
    TV2_s_le = R(s, a_star) + gamma * P_V2

    print("Step 1: Define the difference |T(V1)(s) - T(V2)(s)|.")
    print("Let a_star be the action that maximizes the expression for T(V1)(s).")
    print(f"T(V1)(s) = R(s, a_star) + gamma * {P_V1}")
    print(f"T(V2)(s) <= R(s, a_star) + gamma * sum[P(s'|s,a_star)V2(s')] (because a_star may not be optimal for V2)")
    print("-" * 30)

    print("Step 2: Write the inequality for the difference.")
    print("T(V1)(s) - T(V2)(s) <= (R(s, a_star) + gamma * sum[P(s'|s,a_star)V1(s')]) - (R(s, a_star) + gamma * sum[P(s'|s,a_star)V2(s')])")
    print("-" * 30)
    
    # Symbolic simplification
    difference = TV1_s - TV2_s_le
    
    print("Step 3: Simplify the expression. Notice how the reward term R(s, a_star) cancels out.")
    # In our symbolic representation, this will be explicit.
    print(f"T(V1)(s) - T(V2)(s) <= {difference}")
    simplified_diff = gamma * (P_V1 - P_V2)
    print(f"Which simplifies to: {simplified_diff}")
    print("-" * 30)

    print("Step 4: Take the absolute value and bound the difference.")
    print("|T(V1)(s) - T(V2)(s)| <= gamma * |sum[P(s'|s,a)(V1(s') - V2(s'))]|")
    print("                      <= gamma * sum[P(s'|s,a)|V1(s') - V2(s')|]")
    print("                      <= gamma * sum[P(s'|s,a) * max_s'|V1(s') - V2(s')|]")
    print("                      <= gamma * max_s'|V1(s') - V2(s')| * sum[P(s'|s,a)]")
    print("                      <= gamma * ||V1 - V2||_inf")
    print("-" * 30)

    print("Conclusion: The proof that the operator is a contraction holds regardless of the reward values.")
    print("Therefore, the geometric convergence property of the value iteration sequence does not depend on the rewards being in any specific bounded range.")
    print("The final equation is ||T(V1) - T(V2)||_inf <= gamma * ||V1 - V2||_inf")
    print(f"In this equation, the reward R does not appear.")
    print(f"The number '1' implicitly multiplies the norm terms, and the discount factor is gamma.")


if __name__ == '__main__':
    demonstrate_contraction_proof()