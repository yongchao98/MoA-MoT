import numpy as np

def demonstrate_non_convexity():
    """
    Demonstrates that for d=3, the set of achievable objectives is not convex,
    even with the condition x_i >= 0. This implies scalarization can fail.
    """
    print("--- Analysis for d=3 ---")
    print("We construct a counterexample to show that the objective space is not convex for d=3.")
    print("A non-convex objective space means scalarization is not guaranteed to find the entire Pareto front.")

    # 1. Define the vectors x_i >= 0 for d=3
    x1 = np.array([1.0, 0.0, 0.0])
    x2 = np.array([0.0, 1.0, 0.0])
    x3 = np.array([1.0, 1.0, 1.0])
    print("\nLet's choose n=3 vectors in d=3, all with non-negative components:")
    print(f"x1 = {x1}")
    print(f"x2 = {x2}")
    print(f"x3 = {x3}")

    print("\nThe set of achievable objectives to be maximized is Y = {( (x1^T*w)^2, (x2^T*w)^2, (x3^T*w)^2 ) | ||w|| = 1}.")
    print("This corresponds to Y = { (w1^2, w2^2, (w1+w2+w3)^2) | w1^2+w2^2+w3^2 = 1 }.")

    # 2. Pick two points in the set Y
    w_A = np.array([1.0, 0.0, 0.0])
    p_A = np.array([(x1 @ w_A)**2, (x2 @ w_A)**2, (x3 @ w_A)**2])
    print(f"\nPoint A: For w_A = {w_A}, we get the point P_A = {p_A} in Y.")

    w_B = np.array([0.0, 1.0, 0.0])
    p_B = np.array([(x1 @ w_B)**2, (x2 @ w_B)**2, (x3 @ w_B)**2])
    print(f"Point B: For w_B = {w_B}, we get the point P_B = {p_B} in Y.")

    # 3. Compute their midpoint M and check if M is in Y
    M = (p_A + p_B) / 2
    print(f"\nThe midpoint of P_A and P_B is M = {M}.")

    print("\nCan this point M be generated? We need to find a unit vector w = (w1, w2, w3) such that:")
    print(f"w1^2 = {M[0]}")
    print(f"w2^2 = {M[1]}")
    print(f"(w1+w2+w3)^2 = {M[2]}")

    print("\nThe constraint ||w||^2 = 1 means w1^2 + w2^2 + w3^2 = 1.")
    w3_squared = 1.0 - M[0] - M[1]
    print(f"Using the first two conditions: {M[0]} + {M[1]} + w3^2 = 1  =>  w3^2 = {w3_squared}")
    print("So, any w generating the first two components of M must have w3 = 0.")

    print("\nNow we check the third condition with w3 = 0: (w1 + w2)^2 must be equal to 1.0.")
    w1 = np.sqrt(M[0])
    w2 = np.sqrt(M[1])
    # There are four possibilities for the signs of w1 and w2.
    val1 = (w1 + w2)**2
    val2 = (w1 - w2)**2
    print(f"However, from w1^2=0.5 and w2^2=0.5, the value of (w1+w2)^2 can only be:")
    print(f"({w1:.4f} + {w2:.4f})^2 = {val1}")
    print(f"or ({w1:.4f} - {w2:.4f})^2 = {val2}")

    print(f"\nNeither of these possible values ({val1}, {val2}) is equal to the required value of {M[2]}.")
    print(f"Thus, the midpoint M={M} is not in the set Y.")

    print("\nConclusion: The set of achievable objectives is not convex for d=3.")
    print("Since the property holds for d=2, the largest value of d for which it holds is 2.")

demonstrate_non_convexity()