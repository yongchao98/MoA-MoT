import numpy as np

def demonstrate_lasso_equivalence_issue():
    """
    This function demonstrates that the solution sets for the constrained
    and penalized forms of Lasso are not always identical, particularly
    when the ordinary least squares solution is not unique (e.g., p > N).
    """
    print("Investigating the equivalence of two Lasso formulations.")
    print("Formulation 1 (Constrained): argmin RSS(b) s.t. ||b||_1 <= t")
    print("Formulation 2 (Penalized): argmin RSS(b) + lambda * ||b||_1")
    print("\nWe use an example where p > N to ensure non-unique solutions:")
    print("N=1 observation, p=2 predictors.")
    # y = 1, X = [[1, 1]]
    # The residual sum of squares (RSS) is (1 - (b1*x1 + b2*x2))^2 = (1 - b1 - b2)^2
    print("y = [1], X = [[1, 1]]")
    print("RSS(b) = (1 - b1 - b2)^2")

    # We choose lambda = 1.0. The corresponding t can be shown to be 0.5.
    lam = 1.0
    t = 0.5
    print(f"\nWe will compare the solution sets for lambda = {lam} and t = {t}.")
    print("-" * 50)

    # --- Analysis of the Penalized Problem ---
    # The penalized problem favors sparse solutions due to the nature of the L1 norm.
    # For a fixed sum s = b1+b2, the term |b1|+|b2| is minimized when one coefficient is zero.
    # The optimal sum can be found to be s_opt = 1 - lambda/2.
    s_opt_penalized = 1 - lam / 2
    
    print("1. Penalized Problem Solution Set (lambda = 1.0):")
    print("The problem favors sparse solutions where b1=0 or b2=0.")
    print(f"The optimal sum b1+b2 is {s_opt_penalized}.")
    print("The solution set consists of two points:")
    sol1_penalized = (s_opt_penalized, 0)
    sol2_penalized = (0, s_opt_penalized)
    print(f"  - Solution 1: b = {sol1_penalized}")
    print(f"  - Solution 2: b = {sol2_penalized}")
    print("-" * 50)

    # --- Analysis of the Constrained Problem ---
    # We want to make b1+b2 as close to 1 as possible, s.t. |b1|+|b2| <= t.
    # The maximum value of b1+b2 on the feasible set is t.
    s_opt_constrained = t

    print("2. Constrained Problem Solution Set (t = 0.5):")
    print(f"The optimal sum b1+b2 is {s_opt_constrained}.")
    print("The solution set contains all points (b1, b2) on the line segment connecting (0.5, 0) and (0, 0.5).")
    print("Example points in the solution set:")
    print(f"  - ({t}, 0.0)")
    print(f"  - (0.0, {t})")
    print(f"  - ({t/2}, {t/2})")
    print("-" * 50)

    # --- Conclusion ---
    print("\nConclusion:")
    print("The solution set for the penalized problem contains only two points.")
    print("The solution set for the constrained problem is an entire line segment.")
    print("Since the solution sets are not identical, the equivalence is not strictly true.")
    print("\nThis demonstrates the statement is 'False in general'.")
    print("It is 'often true in practice' because this issue only arises for non-unique solutions, and when solutions are unique, the equivalence holds.")

demonstrate_lasso_equivalence_issue()