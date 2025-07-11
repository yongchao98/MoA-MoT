import itertools

def find_max_matchups():
    """
    Finds the largest list of 4-player matchups from 11 players, where no
    two matchups share more than two players.
    """
    num_players = 11
    group_size = 4
    max_common_players = 2

    players = range(num_players)

    # Generate all possible matchups in a deterministic order
    all_possible_matchups = list(itertools.combinations(players, group_size))

    # This list will hold the final selection of matchups
    selected_matchups = []

    # Implement the greedy algorithm
    for candidate_tuple in all_possible_matchups:
        candidate_matchup = set(candidate_tuple)
        
        is_valid = True
        # Check the candidate against all previously selected matchups
        for existing_matchup in selected_matchups:
            common_players_count = len(candidate_matchup.intersection(existing_matchup))
            
            if common_players_count > max_common_players:
                is_valid = False
                break
        
        # If the candidate is valid with all existing selections, add it
        if is_valid:
            selected_matchups.append(candidate_matchup)

    # Print the results
    print(f"Found {len(selected_matchups)} valid matchups:")
    for matchup in selected_matchups:
        # Sort the players in each matchup for a clean, readable output
        sorted_matchup_list = sorted(list(matchup))
        print(f"Matchup: {sorted_matchup_list[0]}, {sorted_matchup_list[1]}, {sorted_matchup_list[2]}, {sorted_matchup_list[3]}")
    
    print(f"\nThe largest list of matchups found is of size {len(selected_matchups)}.")
    
    # Return the final count for the answer block
    return len(selected_matchups)

# Run the function and print the final answer in the specified format
final_answer = find_max_matchups()
# The final answer is enclosed in <<<>>>
print(f"\n<<<Solution>>>\n<<<The final answer is {final_answer}>>>")
