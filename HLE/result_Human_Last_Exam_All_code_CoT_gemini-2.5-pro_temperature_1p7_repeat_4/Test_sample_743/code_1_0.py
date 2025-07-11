import itertools

# --- Configuration ---
NUM_PLAYERS = 11
GROUP_SIZE = 4
MAX_COMMON_PLAYERS = 2

# Global variable to store the best solution found so far.
# It will be a list of frozensets.
best_solution = []

def is_compatible(candidate_matchup, existing_matchups):
    """
    Checks if a candidate matchup is compatible with all existing matchups in a list.
    A matchup is compatible if it has at most MAX_COMMON_PLAYERS in common.
    """
    for existing in existing_matchups:
        # frozensets allow for efficient intersection
        if len(candidate_matchup.intersection(existing)) > MAX_COMMON_PLAYERS:
            return False
    return True

def find_largest_set_recursive(all_matchups, start_index, current_solution):
    """
    A recursive backtracking function to find the maximum set of compatible matchups.
    
    :param all_matchups: The complete list of C(11, 4) possible matchups.
    :param start_index: The index in all_matchups to start searching from, to avoid duplicates.
    :param current_solution: The list of compatible matchups built so far in this recursive path.
    """
    global best_solution

    # If the current path has yielded a better solution, update the global best.
    if len(current_solution) > len(best_solution):
        best_solution = current_solution

    # Pruning: If the current solution's length plus the remaining number of
    # candidates cannot possibly exceed the best solution found so far,
    # stop exploring this path.
    if len(current_solution) + (len(all_matchups) - start_index) <= len(best_solution):
        return

    # Iterate through the remaining potential matchups
    for i in range(start_index, len(all_matchups)):
        candidate = all_matchups[i]
        
        # Check if this candidate is compatible with all matchups already in our current solution
        if is_compatible(candidate, current_solution):
            # If it is, we add it to the solution for this path and recurse
            find_largest_set_recursive(all_matchups, i + 1, current_solution + [candidate])

def solve_and_print():
    """
    Main function to set up the problem, run the solver, and print the results.
    """
    players = range(NUM_PLAYERS)
    
    # Generate all C(11, 4) = 330 possible matchups.
    # We use frozenset because they are hashable and allow for efficient set operations.
    all_possible_matchups = [frozenset(p) for p in itertools.combinations(players, GROUP_SIZE)]
    
    # Start the recursive search from the beginning.
    find_largest_set_recursive(all_possible_matchups, 0, [])

    # Print the final results. The cryptic instruction "output each number in the final equation!"
    # is interpreted as clearly printing the final result: the size of the set and its elements.
    print(f"The largest list of matchups has a size of {len(best_solution)}.")
    print("One such set of matchups is:")
    
    # Sort the players within each matchup for consistent, readable output.
    for matchup in best_solution:
        print(sorted(list(matchup)))

if __name__ == '__main__':
    solve_and_print()
    print(f'<<<{len(best_solution)}>>>')
