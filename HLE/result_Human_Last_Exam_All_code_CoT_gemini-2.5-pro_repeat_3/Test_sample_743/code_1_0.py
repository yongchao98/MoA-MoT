import itertools

def find_max_matchups():
    """
    This script finds the largest possible list of tennis matchups (groups of 4)
    from 11 players, with the constraint that no two matchups share more than
    two players. It uses a greedy algorithm to build the list of matchups.
    """
    num_players = 11
    group_size = 4
    max_common_players = 2

    # We can represent players with numbers from 1 to 11.
    players = range(1, num_players + 1)

    # 1. Generate all possible matchups.
    all_possible_matchups = list(itertools.combinations(players, group_size))

    # 2. Initialize the list for the final selection.
    final_matchups = []

    # 3. Iterate through all candidates and add them if they are valid.
    for candidate_matchup in all_possible_matchups:
        is_valid = True
        candidate_set = set(candidate_matchup)

        # Check the candidate against all matchups already selected.
        for existing_matchup in final_matchups:
            existing_set = set(existing_matchup)
            
            # Calculate how many players are common between the two matchups.
            common_players_count = len(candidate_set.intersection(existing_set))

            if common_players_count > max_common_players:
                is_valid = False
                break  # This candidate is invalid, so we stop checking.

        # If the candidate is valid with all existing matchups, add it to our list.
        if is_valid:
            final_matchups.append(candidate_matchup)

    # 4. Print the final results.
    print(f"The largest list of matchups found has {len(final_matchups)} entries.\n")
    print("The matchups are:")
    for i, matchup in enumerate(final_matchups):
        # The problem asks to output each number in the final equation.
        # We will print each matchup as a list of player numbers.
        player1, player2, player3, player4 = matchup
        print(f"Matchup {i+1:2d}: Player {player1}, Player {player2}, Player {player3}, Player {player4}")

if __name__ == '__main__':
    find_max_matchups()
