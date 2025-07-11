def solve_utility_problem():
    """
    Analyzes the Three Utilities Problem using Euler's formula for planar graphs.
    """
    # V = Number of vertices
    # There are 3 houses and 3 utilities.
    num_houses = 3
    num_utilities = 3
    V = num_houses + num_utilities
    print(f"Step 1: Calculate the number of vertices (V).")
    print(f"There are {num_houses} houses and {num_utilities} utilities.")
    print(f"V = {num_houses} + {num_utilities} = {V}\n")

    # E = Number of edges
    # Each house connects to each utility.
    E = num_houses * num_utilities
    print(f"Step 2: Calculate the number of edges (E).")
    print(f"Each of the {num_houses} houses connects to all {num_utilities} utilities.")
    print(f"E = {num_houses} * {num_utilities} = {E}\n")

    # F = Number of faces
    # For a planar graph, Euler's formula is V - E + F = 2.
    # We can rearrange it to F = 2 - V + E.
    F = 2 - V + E
    print(f"Step 3: Assume the graph is planar and calculate the required number of faces (F) using Euler's formula (V - E + F = 2).")
    print(f"F = 2 - V + E")
    print(f"F = 2 - {V} + {E} = {F}\n")

    print(f"Step 4: Analyze the properties of the faces.")
    print(f"The graph is bipartite, so the shortest possible cycle has a length of 4.")
    print(f"This means every face must be bounded by at least 4 edges.\n")

    print(f"Step 5: Formulate an inequality based on the face properties.")
    print(f"If F faces each have at least 4 edges, the sum of edges for all faces is at least 4 * F.")
    print(f"Since each edge is shared by two faces, this sum is also equal to 2 * E.")
    print(f"This gives the inequality: 2 * E >= 4 * F, which simplifies to E >= 2 * F.\n")

    print(f"Step 6: Check for a contradiction.")
    print(f"We test our values E = {E} and F = {F} in the inequality E >= 2 * F.")
    print(f"Is {E} >= 2 * {F}?")
    print(f"Is {E} >= {2 * F}?")

    if E >= 2 * F:
        print(f"The condition is met. This does not lead to a contradiction.\n")
    else:
        print(f"The condition is FALSE. {E} is not greater than or equal to {2 * F}.\n")
        print(f"Conclusion: The initial assumption that the graph is planar leads to a mathematical contradiction.")
        print(f"Therefore, the problem is impossible to solve under the given constraints.")
        print(f"This is formally stated by Kuratowski's Theorem, which identifies K(3,3) as a non-planar graph.")

solve_utility_problem()
<<<E>>>