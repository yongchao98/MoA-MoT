import itertools

def find_largest_matchup_list():
    """
    This function solves the problem of finding the largest set of doubles matchups
    from 11 players with the given constraint.

    Problem details:
    - Total players: 11
    - Group size for a matchup: 4
    - Constraint: Any two matchups can share at most 2 players.
    """

    num_players = 11
    group_size = 4
    max_common_players = 2

    # Represent players as numbers 0 through 10
    players = list(range(num_players))

    # Step 1: Generate all possible unique matchups (combinations)
    all_possible_matchups = list(itertools.combinations(players, group_size))
    
    # Step 2: Use a greedy approach to build the list of valid matchups
    solution_matchups = []

    # Iterate through each potential matchup
    for candidate_matchup in all_possible_matchups:
        candidate_set = set(candidate_matchup)
        is_valid = True
        
        # Check the candidate against all matchups already selected for our solution
        for existing_matchup in solution_matchups:
            existing_set = set(existing_matchup)
            
            # Calculate how many players are common to both matchups
            common_players_count = len(candidate_set.intersection(existing_set))
            
            # If the number of common players exceeds the allowed maximum,
            # this candidate is not valid.
            if common_players_count > max_common_players:
                is_valid = False
                break  # Exit the inner loop; no need to check further
        
        # If the candidate is valid against all existing solution matchups, add it.
        if is_valid:
            solution_matchups.append(candidate_matchup)
            
    # Step 3: Print the final results
    print(f"There are {num_players} players to be arranged in groups of {group_size}.")
    print(f"The constraint is that no two groups can have more than {max_common_players} players in common.")
    print("-" * 50)
    print(f"The largest list of matchups that can be created is: {len(solution_matchups)}")
    print("\nThe list of matchups is (players are numbered 0 to 10):")
    
    # Output each matchup in the final list
    for i, matchup in enumerate(solution_matchups):
        # We add 1 to each player number for a more human-readable "Player 1" to "Player 11" format,
        # but the actual numbers in the calculation are 0-10.
        # Let's stick to the 0-10 format for consistency with the algorithm's logic.
        print(f"Matchup {i+1}: {list(matchup)}")

# Execute the function to solve the problem
find_largest_matchup_list()
<<<12>>>