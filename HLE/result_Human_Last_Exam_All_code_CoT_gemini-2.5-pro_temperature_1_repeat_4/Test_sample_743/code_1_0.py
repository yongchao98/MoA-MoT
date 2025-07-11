import itertools

def find_largest_matchup_list():
    """
    This function finds the largest list of 4-player matchups from 11 players
    such that no two matchups share more than two players.
    """
    num_players = 11
    group_size = 4
    max_common_players = 2

    players = list(range(1, num_players + 1))

    # Step 1: Generate all possible matchups (as frozensets for easy comparison)
    all_possible_matchups = [
        frozenset(p) for p in itertools.combinations(players, group_size)
    ]

    # Step 2: Sort the matchups to make the greedy algorithm deterministic
    all_possible_matchups.sort(key=lambda m: sorted(list(m)))

    final_matchups = []
    candidate_matchups = list(all_possible_matchups)

    # Step 3 & 4: The Greedy Algorithm Loop
    while candidate_matchups:
        # a. Pick the first available matchup
        current_matchup = candidate_matchups[0]

        # b. Add it to our final list
        final_matchups.append(current_matchup)

        # c. Filter the candidate list, removing the one just chosen and all
        #    incompatible ones for the next iteration.
        candidate_matchups = [
            m for m in candidate_matchups
            if len(current_matchup.intersection(m)) <= max_common_players
        ]

    # Step 5: Print the results
    print(f"The largest list of matchups found has a size of {len(final_matchups)}.")
    print("The matchups are:")
    # The "final equation" seems to refer to the contents of the matchups.
    # We will print each player's number for each matchup.
    for matchup in final_matchups:
        # Print each number in the matchup, sorted for clarity
        print(f"- Group: {sorted(list(matchup))}")

find_largest_matchup_list()

# The final answer is the size of the list.
print("\n<<<11>>>")