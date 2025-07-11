import itertools
import networkx as nx

def solve_tennis_matchups():
    """
    Finds the largest list of tennis matchups satisfying the given condition
    by modeling the problem as finding the maximum clique in a graph.
    """
    players = range(11)
    
    # Step 1 & 2: Generate all 330 possible groups of 4 from 11 players.
    # These are our potential matchups and will become the nodes of our graph.
    all_possible_matchups = list(itertools.combinations(players, 4))
    all_matchups_sets = [frozenset(m) for m in all_possible_matchups]
    num_nodes = len(all_matchups_sets)

    # Step 3 & 4: Create the "compatibility graph".
    # An edge exists between two matchups if they share 2 or fewer players.
    G = nx.Graph()
    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):
            intersection_size = len(all_matchups_sets[i].intersection(all_matchups_sets[j]))
            if intersection_size <= 2:
                G.add_edge(i, j)

    # Step 5 & 6: Find the maximum clique. This corresponds to the largest set
    # of mutually compatible matchups.
    # Note: This computation might take a few seconds.
    max_clique_indices = max(nx.find_cliques(G), key=len)

    largest_size = len(max_clique_indices)
    
    print(f"The largest list of matchups that can be created is of size: {largest_size}")
    
    # As requested, outputting the matchups found. Each matchup is a collection
    # of numbers representing the players.
    print("\nOne such list of matchups (represented by player indices 0-10) is:")
    final_matchups = sorted([all_possible_matchups[i] for i in max_clique_indices])

    for i, matchup in enumerate(final_matchups):
        # The players in each matchup are printed below.
        print(f"Matchup {i+1}: {matchup}")

# Execute the plan.
solve_tennis_matchups()