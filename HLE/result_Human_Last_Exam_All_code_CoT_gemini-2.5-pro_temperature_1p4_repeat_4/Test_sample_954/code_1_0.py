import math

def solve_and_explain():
    """
    This function explains the reasoning step-by-step and prints the final answer.
    """

    print("### Analysis of the Graph Process ###")

    print("\nStep 1: Interpreting the Process Dynamics")
    print("The process describes life reduction for vertices. In each step 'i', for every pair of adjacent and alive vertices (u, v), they both lose min(1/d_u^i, 1/d_v^i) life. The total life lost by a vertex is the sum of losses over all its neighbors.")
    print("Life Lost(u) at step i = sum_{v ~ u} min(1/d_u^i, 1/d_v^i)")

    print("\nStep 2: Identifying Key Vertex Removals")
    print("Consider a vertex 'u' that is a local maximum in terms of degree, i.e., d_u >= d_v for all its neighbors 'v'.")
    print("For such a vertex, min(1/d_u, 1/d_v) = 1/d_u for all neighbors.")
    print("The life it loses is: sum_{v ~ u} (1/d_u) = d_u * (1/d_u) = 1.")
    print("Since each vertex starts with 1 life, any such vertex is removed in a single step where it satisfies this condition.")

    print("\nStep 3: Bounding the Number of Steps (T)")
    print("At any step 'k', let Delta_k be the maximum degree among all currently alive vertices.")
    print("Any vertex 'u' with degree Delta_k is a local maximum degree vertex because d_u = Delta_k >= d_v for all its neighbors 'v'.")
    print("Therefore, all vertices with the current maximum degree are removed in that step.")
    print("This means the maximum degree of the graph of surviving vertices strictly decreases in every step: Delta_{k+1} < Delta_k.")
    print("This implies the total number of steps T is at most the initial maximum degree: T <= Delta.")
    print("Furthermore, it is possible to construct a 'caterpillar' graph where the number of steps is Omega(Delta).")
    print("Thus, the maximum number of steps for a graph with maximum degree Delta is T = Theta(Delta).")

    print("\nStep 4: Analyzing the Three Cases")

    # Case 1
    print("\n--- Case 1: Maximum degree Delta <= sqrt(log n) ---")
    print("The maximum number of steps is f_1(n) = Theta(Delta) = Theta(sqrt(log n)).")
    print("Let's compare this to the given categories. Let x = log n and y = log(log n).")
    print("f_1(n) = Theta(sqrt(x)) = Theta(exp(0.5*y)).")
    print("Category 5 is 2^O(sqrt(log log n)) = exp(O(sqrt(y))). f_1(n) grows faster than this.")
    print("Category 6 is 2^omega(sqrt(log log n)) and O((log n)^0.9). f_1(n) = exp(0.5y) is exp(omega(sqrt(y))).")
    print("And sqrt(log n) is O((log n)^0.9). So, f_1(n) fits in Category 6.")
    f1 = 6

    # Case 2
    print("\n--- Case 2: Maximum degree Delta <= log n ---")
    print("The maximum number of steps is f_2(n) = Theta(Delta) = Theta(log n).")
    print("This directly matches Category 8: f(n) = Theta(log n).")
    f2 = 8

    # Case 3
    print("\n--- Case 3: Any forest ---")
    print("For any forest on n vertices, we want to find the maximum possible number of steps.")
    print("T(G) = Theta(Delta(G)). We need to find the graph structure that maximizes T.")
    print("The caterpillar-like construction that gives T = Omega(Delta) requires n = O(Delta^2) vertices. This means the maximum Delta achievable for this high-step-count structure on n vertices is Delta = Theta(sqrt(n)).")
    print("This leads to a maximum number of steps f_3(n) = Theta(sqrt(n)).")
    print("Let's analyze f_3(n) = Theta(sqrt(n)) = Theta(exp(0.5 * log n)).")
    print("This function grows faster than any polynomial of log n, so it is omega(log n).")
    print("This fits into Category 9: f(n) = omega(log n).")
    f3 = 9

    print("\n### Final Result ###")
    print(f"The digit for f_1(n) is: {f1}")
    print(f"The digit for f_2(n) is: {f2}")
    print(f"The digit for f_3(n) is: {f3}")
    final_number = f"{f1}{f2}{f3}"
    print(f"\nThe resulting three-digit number is {final_number}.")

# Run the analysis
solve_and_explain()