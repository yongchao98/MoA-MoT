import itertools

def find_largest_matchup_list():
    """
    Finds the largest list of 4-player matchups from 11 players such that
    no two matchups have more than two players in common.
    """
    num_players = 11
    group_size = 4
    max_common_players = 2

    players = list(range(1, num_players + 1))

    # Step 1: Generate all possible matchups (combinations of 4 players)
    all_possible_matchups = list(itertools.combinations(players, group_size))

    # Step 2: Use a greedy algorithm to build the list of valid matchups
    final_matchups = []

    for matchup in all_possible_matchups:
        # Assume the current matchup is valid to add
        is_valid_to_add = True

        # Check this matchup against all matchups already added to our final list
        for existing_matchup in final_matchups:
            # Find the number of common players between the two matchups
            common_players = set(matchup).intersection(set(existing_matchup))
            
            # If the number of common players exceeds the allowed maximum, this matchup cannot be added
            if len(common_players) > max_common_players:
                is_valid_to_add = False
                break  # Move to the next potential matchup

        # If the matchup is valid with respect to the entire existing list, add it
        if is_valid_to_add:
            final_matchups.append(matchup)

    # Step 3: Print the final result
    print(f"The largest list of matchups that can be created is of size {len(final_matchups)}.")
    print("The matchups are:")
    for matchup in final_matchups:
        # The instruction "output each number in the final equation" is interpreted
        # as printing each group of players.
        print(f"Group: {list(matchup)}")

# Execute the function to find and print the solution
find_largest_matchup_list()