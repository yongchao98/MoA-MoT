import itertools

def solve_tennis_matchups():
    """
    Finds the largest list of 4-player matchups from 11 players such that
    no two matchups share more than two players.
    """
    num_players = 11
    group_size = 4
    max_common_players = 2

    players = list(range(1, num_players + 1))

    # Generate all possible matchups (combinations of 4 from 11)
    all_possible_matchups = list(itertools.combinations(players, group_size))

    # Use a greedy approach to build the list of valid matchups
    final_matchups = []

    for candidate_matchup in all_possible_matchups:
        is_valid = True
        # Convert to set for efficient intersection
        candidate_set = set(candidate_matchup)
        
        for existing_matchup in final_matchups:
            existing_set = set(existing_matchup)
            
            # Find the number of common players
            common_players = len(candidate_set.intersection(existing_set))
            
            if common_players > max_common_players:
                is_valid = False
                break  # This candidate conflicts with an existing matchup

        if is_valid:
            final_matchups.append(candidate_matchup)

    # Print the final result
    print(f"The largest list contains {len(final_matchups)} matchups.")
    print("The matchups are:")
    for i, matchup in enumerate(final_matchups, 1):
        # Formatting to print each player in the matchup
        player_numbers = " ".join(map(str, matchup))
        print(f"Match {i:2d}: [ {player_numbers} ]")

if __name__ == '__main__':
    solve_tennis_matchups()
