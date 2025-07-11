import itertools

def solve_tennis_matchups():
    """
    Solves the tennis matchup problem by finding the maximum number of
    groups of 4 from 11 players, where no two groups share more than 2 players.

    This is equivalent to finding the maximum independent set on a graph
    where vertices are matchups and edges connect conflicting matchups.
    """
    NUM_PLAYERS = 11
    GROUP_SIZE = 4
    MAX_COMMON_PLAYERS = 2

    players = range(NUM_PLAYERS)
    
    # 1. Generate all possible matchups
    all_matchups = list(itertools.combinations(players, GROUP_SIZE))
    num_matchups = len(all_matchups)

    # 2. Build the conflict graph
    # conflicts[i] is a set of indices j > i where matchup i and j conflict.
    conflicts = [set() for _ in range(num_matchups)]
    for i in range(num_matchups):
        for j in range(i + 1, num_matchups):
            # A conflict occurs if the intersection size is greater than MAX_COMMON_PLAYERS
            intersection_size = len(set(all_matchups[i]) & set(all_matchups[j]))
            if intersection_size > MAX_COMMON_PLAYERS:
                conflicts[i].add(j)

    # 3. Backtracking search for the maximum independent set
    max_solution = []
    
    # The recursive function to explore solutions
    def find_max(candidate_indices, current_solution):
        nonlocal max_solution
        
        # Pruning: if the current path can't possibly beat the best found so far, stop.
        if len(current_solution) + len(candidate_indices) <= len(max_solution):
            return

        # Iterate through the candidates to extend the current solution
        for i, idx in enumerate(candidate_indices):
            
            # Check if the new candidate conflicts with the current solution
            is_compatible = True
            for sol_idx in current_solution:
                # Use the pre-computed conflict graph for efficiency
                if sol_idx < idx and idx in conflicts[sol_idx]:
                    is_compatible = False
                    break
                # The other case is not needed due to how the graph is built (j > i)
            
            if is_compatible:
                # If compatible, recurse with the new candidate added to the solution
                new_candidates = [c for c in candidate_indices[i+1:] if c not in conflicts[idx]]
                find_max(new_candidates, current_solution + [idx])

        # After exploring all extensions from the current state,
        # check if the current solution itself is a new maximum
        if len(current_solution) > len(max_solution):
            max_solution = current_solution
            print(f"Found a new largest list with {len(max_solution)} matchups.")

    # Start the search with all matchups as initial candidates
    find_max(list(range(num_matchups)), [])

    # Print the final result
    print("\n----------------------------------------------------")
    print("The largest possible list of matchups has been found.")
    print(f"Maximum number of matchups: {len(max_solution)}")
    print("The matchups are:")
    for i, matchup_idx in enumerate(max_solution):
        # Format players to be 1-indexed for readability if preferred, but 0-indexed is fine.
        players_in_matchup = all_matchups[matchup_idx]
        print(f"  Matchup {i+1}: Players {players_in_matchup}")

solve_tennis_matchups()