import itertools

def find_largest_matchup_list():
    """
    This function solves the tennis matchup problem using a backtracking algorithm.
    It finds the largest set of matchups (groups of 4 from 11 players)
    where any two matchups share at most two players.
    """
    players = range(1, 12)  # We have 11 players, numbered 1 to 11
    group_size = 4
    max_shared_players = 2

    # 1. Generate all possible matchups as sets for easy intersection.
    all_matchups = [set(c) for c in itertools.combinations(players, group_size)]
    
    # This list will store the largest valid set of matchups found.
    best_solution = []

    def backtrack(start_index, current_solution):
        """
        The recursive backtracking function.
        - start_index: The index in all_matchups to start searching from.
        - current_solution: The list of valid matchups found so far in this path.
        """
        nonlocal best_solution

        # --- Optimization ---
        # Prune the search if the current path can't possibly beat the best solution found so far.
        if len(current_solution) + (len(all_matchups) - start_index) <= len(best_solution):
            return

        # If the current solution is better than our best so far, update the best solution.
        if len(current_solution) > len(best_solution):
            best_solution = list(current_solution)
        
        # Iterate through the remaining matchups to try and extend our solution.
        for i in range(start_index, len(all_matchups)):
            new_matchup = all_matchups[i]
            
            # Check if this new matchup is compatible with our current solution.
            is_compatible = True
            for existing_matchup in current_solution:
                if len(new_matchup.intersection(existing_matchup)) > max_shared_players:
                    is_compatible = False
                    break
            
            # If it's compatible, add it and continue the search from the next index.
            if is_compatible:
                current_solution.append(new_matchup)
                backtrack(i + 1, current_solution)
                current_solution.pop()  # Backtrack: remove the matchup to explore other paths.

    # Start the search from the beginning.
    backtrack(0, [])

    # Print the final results.
    print(f"The largest list of matchups that can be created is: {len(best_solution)}")
    print("\nThe matchups are:")
    for matchup in best_solution:
        # Sort for consistent output and format nicely
        players_in_matchup = sorted(list(matchup))
        print(", ".join(map(str, players_in_matchup)))

if __name__ == '__main__':
    find_largest_matchup_list()