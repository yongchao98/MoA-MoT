import itertools

def solve_tennis_matchups():
    """
    Solves the tennis matchup problem by modeling it as a maximum clique problem.
    """
    
    # 1. Define problem parameters
    NUM_PLAYERS = 11
    GROUP_SIZE = 4
    MAX_COMMON_PLAYERS = 2

    # 2. Generate all possible matchups (vertices of the graph)
    players = range(NUM_PLAYERS)
    all_matchups = list(itertools.combinations(players, GROUP_SIZE))

    print("Step 1: Generated all {} possible matchups.".format(len(all_matchups)))

    # 3. Build the compatibility graph
    # An edge exists between two matchups if they are compatible.
    graph = {matchup: set() for matchup in all_matchups}
    for i in range(len(all_matchups)):
        for j in range(i + 1, len(all_matchups)):
            m1 = all_matchups[i]
            m2 = all_matchups[j]
            
            common_players = len(set(m1).intersection(set(m2)))
            
            if common_players <= MAX_COMMON_PLAYERS:
                graph[m1].add(m2)
                graph[m2].add(m1)

    print("Step 2: Built the compatibility graph.")
    
    # This list will store the largest clique found.
    max_clique_found = []

    def find_max_clique(R, P, X):
        """
        Bron-Kerbosch algorithm with pivoting to find the maximum clique.
        R: The set of vertices in the current clique.
        P: The set of candidate vertices to extend the clique.
        X: The set of vertices already processed.
        """
        nonlocal max_clique_found
        
        if not P and not X:
            # R is a maximal clique. Check if it's the largest.
            if len(R) > len(max_clique_found):
                max_clique_found = list(R)
            return

        if not P:
            return

        # Choose a pivot vertex 'u' from P union X to optimize.
        # We select the vertex with the most neighbors in P.
        try:
            pivot_candidates = P | X
            pivot = max(pivot_candidates, key=lambda u: len(P.intersection(graph[u])))
        except (ValueError, KeyError):
            return

        # Iterate through candidates not connected to the pivot
        P_without_neighbors_of_pivot = P - graph[pivot]

        for v in list(P_without_neighbors_of_pivot):
            find_max_clique(
                R | {v},
                P.intersection(graph[v]),
                X.intersection(graph[v])
            )
            P.remove(v)
            X.add(v)

    print("Step 3: Searching for the largest list of compatible matchups...")
    
    # 4. Run the algorithm on the graph
    # Initial call with R and X empty, and P containing all vertices.
    find_max_clique(set(), set(all_matchups), set())
    
    # 5. Print the final result
    print("\n--- Solution ---")
    final_size = len(max_clique_found)
    print("The largest list of matchups that can be created is: {}".format(final_size))
    print("\nOne possible list of these matchups is:")
    
    # Sort for consistent output format
    sorted_matchups = sorted([list(matchup) for matchup in max_clique_found])

    for matchup in sorted_matchups:
        # The problem asks to output each number in the final list
        print("Matchup: {}".format(matchup))

solve_tennis_matchups()
<<<11>>>