def check_graph_properties(n):
    """
    Analyzes the mathematical constraints on a graph with n vertices based on the problem statement.

    The properties are:
    1. n vertices
    2. 7-regular
    3. Chromatic number is 5
    4. Exactly n copies of C5
    5. No three C5s share a common vertex

    This function demonstrates that these properties are self-contradictory.
    """
    print(f"Analyzing for n = {n}")
    print("-" * 20)
    print("Let n_k be the number of vertices belonging to exactly k C5-cycles.")
    print("From the problem constraints, a vertex can belong to 0, 1, or 2 cycles.")
    print("So, the total number of vertices is n = n_0 + n_1 + n_2.")
    print("\nLet's count the total (vertex, C5) incidences in two ways:")
    print(f"1. Per cycle: n cycles * 5 vertices/cycle = {5*n} incidences.")
    print("2. Per vertex: (n_0 * 0) + (n_1 * 1) + (n_2 * 2) incidences.")
    print("This gives us two equations:")
    print("Eq 1: n_0 + n_1 + n_2 = n")
    print("Eq 2: n_1 + 2*n_2 = 5*n")

    print("\nSolving these equations:")
    print("From Eq 1, we get n_1 = n - n_0 - n_2.")
    print("Substituting n_1 into Eq 2: 5*n = (n - n_0 - n_2) + 2*n_2")
    print("Simplifying: 5*n = n - n_0 + n_2")
    print("This gives a formula for n_2: n_2 = 4*n + n_0")

    print("\nChecking for contradictions:")
    print(f"For n = {n}, n_2 must be equal to 4*{n} + n_0, which is {4*n} + n_0.")
    print("However, n_2 is the number of vertices of a certain type, so n_2 cannot be greater than the total number of vertices, n.")
    print(f"So, we must have n_2 <= n, which means {4*n} + n_0 <= {n}.")
    print(f"This simplifies to the inequality: {3*n} + n_0 <= 0.")

    three_n = 3 * n
    print("\nFinal conclusion:")
    print(f"For n={n}, 3*n = {three_n}. Since n must be a positive integer, 3*n is positive.")
    print("n_0 (number of vertices) must be non-negative (>= 0).")
    print(f"The expression '{three_n} + n_0' is therefore the sum of a positive and a non-negative number, which must be positive.")
    print(f"The condition '{three_n} + n_0 <= 0' can never be true.")
    print("This proves that no such graph exists for this n.")


# The problem asks for the smallest composite n. Smallest composite is 4.
# A 7-regular graph needs at least 8 vertices, and an odd-degree regular graph must have an even number of vertices.
# Let's test the smallest possible composite candidate, n=8.
smallest_possible_n = 8
check_graph_properties(smallest_possible_n)
