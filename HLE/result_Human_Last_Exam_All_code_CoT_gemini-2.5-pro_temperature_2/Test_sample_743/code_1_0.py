import itertools

# Configuration for the problem
NUM_PLAYERS = 11
GROUP_SIZE = 4

def solve_tennis_matchups():
    """
    Finds the largest list of tennis matchups (groups of 4 from 11 players)
    such that no two matchups share more than two players.
    """
    print("Thinking...")
    print(f"Generating all possible matchups of size {GROUP_SIZE} from {NUM_PLAYERS} players...")
    
    players = range(NUM_PLAYERS)
    
    # Using frozenset for matchups because they are hashable and can be members of sets.
    all_matchups = list(itertools.combinations(players, GROUP_SIZE))
    all_matchups = [frozenset(m) for m in all_matchups]
    num_matchups = len(all_matchups)
    
    print(f"Total possible matchups: {num_matchups}")
    print("Building the compatibility matrix... (this may take a moment)")

    # Pre-compute which pairs of matchups are incompatible
    # is_incompatible[i][j] is True if matchup i and j have an intersection > 2
    is_incompatible = [[False] * num_matchups for _ in range(num_matchups)]
    for i in range(num_matchups):
        for j in range(i + 1, num_matchups):
            m1 = all_matchups[i]
            m2 = all_matchups[j]
            if len(m1.intersection(m2)) > 2:
                is_incompatible[i][j] = True
                is_incompatible[j][i] = True

    # Global variable to store the indices of the matchups in the best list found
    best_solution_indices = []

    def find_max_clique(potential_indices, current_clique_indices):
        """
        A recursive backtracking function to find the maximum clique.
        
        Args:
            potential_indices: Indices of matchups that are compatible with the current clique
                               and can be added.
            current_clique_indices: The list of indices forming the current clique.
        """
        nonlocal best_solution_indices

        # Pruning optimization: if the current clique plus all remaining candidates
        # can't beat the best solution found so far, we stop exploring this path.
        if len(current_clique_indices) + len(potential_indices) <= len(best_solution_indices):
            return

        # When we run out of potential candidates, the current clique is maximal.
        # We check if it's the new best overall solution.
        if not potential_indices:
            if len(current_clique_indices) > len(best_solution_indices):
                best_solution_indices = current_clique_indices[:]
            return
            
        # Iterate through the potential candidates to extend the current clique
        for i in range(len(potential_indices)):
            pivot_candidate_index = potential_indices[i]
            
            # Form a new list of candidates for the recursive call.
            # These must be compatible with the pivot candidate we're adding.
            new_potential_indices = []
            for j in range(i + 1, len(potential_indices)):
                next_potential_index = potential_indices[j]
                if not is_incompatible[pivot_candidate_index][next_potential_index]:
                    new_potential_indices.append(next_potential_index)
            
            # Recursive call
            find_max_clique(new_potential_indices, current_clique_indices + [pivot_candidate_index])

    print("Searching for the largest list of compatible matchups...")
    
    # Start the search with all matchups as initial candidates
    find_max_clique(list(range(num_matchups)), [])
    
    # Output the final results
    print("\n--- Search Complete ---")
    
    max_size = len(best_solution_indices)
    print(f"\nThe largest list of matchups that can be created has {max_size} matchups.")
    print("Here is one such list:")
    
    # Get the actual matchups from their indices
    final_matchup_list = [all_matchups[i] for i in best_solution_indices]
    
    for i, matchup in enumerate(final_matchup_list, 1):
        # Sort players in the matchup for consistent display
        sorted_players = tuple(sorted(list(matchup)))
        print(f"Matchup {i:>2}: {sorted_players}")

    print(f"\n<<<{max_size}>>>")

if __name__ == '__main__':
    solve_tennis_matchups()