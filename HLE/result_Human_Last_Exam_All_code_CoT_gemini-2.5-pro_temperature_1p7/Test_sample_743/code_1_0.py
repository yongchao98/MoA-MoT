import itertools

def find_largest_matchup_list():
    """
    Finds the largest list of 4-player matchups from 11 players, where no two
    matchups share more than two players.
    """
    num_players = 11
    group_size = 4
    max_common_players = 2

    players = list(range(1, num_players + 1))

    # Step 1: Generate all possible matchups
    all_possible_matchups = list(itertools.combinations(players, group_size))

    # Step 2: Start with an empty list for the final selection of matchups
    final_matchups = []

    # Step 3 & 4: Iterate through all matchups and add them if they are valid
    for new_matchup in all_possible_matchups:
        is_valid_to_add = True
        
        # Check the new matchup against all previously selected matchups
        for existing_matchup in final_matchups:
            # Find how many players are common between the two matchups
            common_players = set(new_matchup).intersection(set(existing_matchup))
            
            # If the number of common players exceeds the limit, this new matchup cannot be added
            if len(common_players) > max_common_players:
                is_valid_to_add = False
                break  # Move on to the next potential matchup
        
        # If the new matchup is valid with respect to all existing ones, add it to our list
        if is_valid_to_add:
            final_matchups.append(new_matchup)

    # Step 5: Print the final results
    print(f"The largest list of matchups that can be created is: {len(final_matchups)}")
    print("\nThe matchups in this list are:")
    for matchup in final_matchups:
        # Convert tuple to list for printing
        print(list(matchup))

# Run the function to find and print the solution
find_largest_matchup_list()