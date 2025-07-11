import itertools

def find_max_matchups():
    """
    This function calculates the largest list of 4-player matchups from 11 players
    such that no two matchups share more than two players.
    """
    num_players = 11
    matchup_size = 4
    players = range(num_players)

    # Generate all possible matchups (combinations of 4 from 11)
    all_possible_matchups = list(itertools.combinations(players, matchup_size))

    # This is a greedy algorithm to build the list of final matchups
    final_matchups = []

    for candidate_matchup in all_possible_matchups:
        is_compatible = True
        candidate_set = set(candidate_matchup)

        # Check if the candidate is compatible with all previously chosen matchups
        for final_matchup in final_matchups:
            final_set = set(final_matchup)
            
            # Check the number of common players (intersection size)
            if len(candidate_set.intersection(final_set)) > 2:
                is_compatible = False
                break  # Conflict found, stop checking this candidate

        # If the candidate is compatible with all existing final matchups, add it
        if is_compatible:
            final_matchups.append(candidate_matchup)

    # Print the final result
    print(f"The largest list of matchups that can be created is: {len(final_matchups)}")
    print("One such list of matchups is:")
    for matchup in final_matchups:
        # The prompt asks to output each number in the final equation.
        # We interpret this as listing the players in each valid matchup.
        print(f"Matchup with players: {matchup[0]}, {matchup[1]}, {matchup[2]}, {matchup[3]}")

if __name__ == '__main__':
    find_max_matchups()