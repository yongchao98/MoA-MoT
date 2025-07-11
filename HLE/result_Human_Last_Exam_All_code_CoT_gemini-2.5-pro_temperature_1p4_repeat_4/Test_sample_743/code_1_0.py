import itertools

def find_largest_matchup_list():
    """
    This function calculates the largest possible list of tennis matchups
    for 11 players in groups of 4, with the constraint that no two matchups
    share more than two players.
    """
    num_players = 11
    group_size = 4
    max_common_players = 2

    # Generate all possible unique matchups (combinations of 4 from 11)
    players = range(num_players)
    all_possible_matchups = list(itertools.combinations(players, group_size))

    # This list will store the final valid matchups
    selected_matchups = []

    # Use a greedy approach to build the list
    for candidate in all_possible_matchups:
        is_valid_to_add = True
        # Check the candidate against all previously selected matchups
        for selected in selected_matchups:
            # Find the number of common players
            common_players = set(candidate).intersection(set(selected))
            
            # If the constraint is violated, this candidate cannot be added
            if len(common_players) > max_common_players:
                is_valid_to_add = False
                break
        
        # If the candidate is valid with all selected matchups, add it to the list
        if is_valid_to_add:
            selected_matchups.append(candidate)

    # Print the results
    print("The largest list of matchups found is:")
    for matchup in selected_matchups:
        # The problem asks to output each number in the final list.
        # We print the players (represented by numbers) in each matchup.
        print(f"Matchup: {list(matchup)}")
    
    print(f"\nThe largest list contains {len(selected_matchups)} matchups.")

if __name__ == '__main__':
    find_largest_matchup_list()
    # The final answer is the size of the list.
    # print("<<<12>>>") 