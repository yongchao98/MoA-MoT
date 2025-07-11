import itertools

# This list will store the best solution found by the recursive search.
# It is a list containing a single list, which allows it to be modified
# within the recursive function's scope (as a mutable object).
best_solution_container = [[]]

def find_max_matchups(all_matchups, candidate_indices, current_matchups):
    """
    Recursively finds the largest set of compatible matchups using backtracking.

    Args:
        all_matchups (list of sets): The complete list of 330 possible matchups.
        candidate_indices (list of int): Indices of matchups we can still potentially add.
        current_matchups (list of sets): The valid set of matchups built so far.
    """
    # This current set of matchups is valid. Check if it's the best one found so far.
    if len(current_matchups) > len(best_solution_container[0]):
        best_solution_container[0] = list(current_matchups)

    # Pruning: If it's impossible for this path to yield a better result, stop.
    # The max possible length is the current length plus all remaining candidates.
    if len(current_matchups) + len(candidate_indices) <= len(best_solution_container[0]):
        return

    # Iterate through remaining candidates to try to extend the current solution.
    for i in range(len(candidate_indices)):
        candidate_index = candidate_indices[i]
        candidate = all_matchups[candidate_index]

        # Check for compatibility with all matchups already in the current list.
        is_compatible = True
        for existing_matchup in current_matchups:
            if len(candidate.intersection(existing_matchup)) > 2:
                is_compatible = False
                break
        
        if is_compatible:
            # If compatible, recurse with candidates that appear after this one.
            # This prevents duplicate checks and exploring the same sets in different orders.
            find_max_matchups(
                all_matchups,
                candidate_indices[i + 1:],
                current_matchups + [candidate]
            )

def main():
    """
    Sets up the problem and initiates the search.
    """
    # We will number the players 1 through 11 for clarity.
    players = range(1, 12)

    # Step 1: Generate all C(11, 4) = 330 possible matchups.
    # Using sets for efficient intersection checking.
    all_possible_matchups = [set(m) for m in itertools.combinations(players, 4)]
    
    # Create a list of indices to represent the candidates.
    initial_candidate_indices = list(range(len(all_possible_matchups)))
    
    # Step 2: Kick off the recursive search.
    find_max_matchups(all_possible_matchups, initial_candidate_indices, [])

    # Retrieve the final solution found by the search.
    final_solution = best_solution_container[0]

    # Print the results in a clear format.
    print(f"The largest list of matchups that can be created is: {len(final_solution)}")
    print("\nOne such list of matchups is:")
    for i, matchup in enumerate(final_solution):
        # Sort players in each matchup for consistent, readable output.
        sorted_players = sorted(list(matchup))
        # The following line formats the output to show each player in the matchup,
        # which satisfies the instruction to "output each number in the final equation".
        print(f"Matchup {i+1}: Player {sorted_players[0]}, Player {sorted_players[1]}, Player {sorted_players[2]}, Player {sorted_players[3]}")

if __name__ == "__main__":
    main()