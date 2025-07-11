import itertools

def find_max_matchups():
    """
    Finds the largest list of 4-player matchups from 11 players,
    where any two matchups share at most 2 players.
    """
    num_players = 11
    group_size = 4
    players = range(1, num_players + 1)

    # 1. Generate all possible matchups
    all_possible_matchups = list(itertools.combinations(players, group_size))

    # 2. Start with an empty list for the final valid matchups
    final_matchups = []

    # 3. Iterate through all candidates and build the list greedily
    for candidate_matchup in all_possible_matchups:
        is_valid = True
        # 4. Check the candidate against all matchups already selected
        for existing_matchup in final_matchups:
            # Find the number of common players
            common_players = set(candidate_matchup).intersection(set(existing_matchup))
            if len(common_players) > 2:
                is_valid = False
                break  # This candidate is invalid, no need to check further
        
        # 5. If the candidate is valid, add it to our list
        if is_valid:
            final_matchups.append(candidate_matchup)

    # Print the results
    print("The largest list of matchups found is:")
    for i, matchup in enumerate(final_matchups, 1):
        # The instruction "output each number in the final equation" is interpreted
        # as showing the elements that form the final result.
        # Here we print each matchup that constitutes the final list.
        print(f"Matchup {i}: {matchup[0]}, {matchup[1]}, {matchup[2]}, {matchup[3]}")

    print("\n-------------------------------------------------")
    print(f"The total number of matchups is: {len(final_matchups)}")
    print("-------------------------------------------------")


if __name__ == '__main__':
    find_max_matchups()
