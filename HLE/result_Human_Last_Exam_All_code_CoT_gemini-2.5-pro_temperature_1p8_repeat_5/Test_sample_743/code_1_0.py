import itertools

def solve_tennis_matchups():
    """
    Finds the largest list of tennis matchups (groups of 4) from 11 players,
    such that no two matchups share more than two players.
    """
    players = list(range(1, 12))

    # 1. Generate all possible matchups of 4 players.
    all_possible_matchups = list(itertools.combinations(players, 4))

    # 2. Use a greedy approach to build the list of valid matchups.
    final_matchups = []

    for candidate_matchup in all_possible_matchups:
        is_compatible = True
        
        # 3. Check compatibility with already selected matchups.
        for existing_matchup in final_matchups:
            # Calculate how many players are common between the two matchups.
            # We convert the tuples to sets to use the intersection operation.
            intersection_size = len(set(candidate_matchup).intersection(set(existing_matchup)))
            
            # If they share more than 2 players, they are not compatible.
            if intersection_size > 2:
                is_compatible = False
                break
        
        # 4. If the candidate is compatible, add it to our final list.
        if is_compatible:
            final_matchups.append(candidate_matchup)

    # 5. Print the results.
    print("The largest list of matchups found is:")
    for i, matchup in enumerate(final_matchups):
        # Format the output to be readable, showing each matchup.
        player_list_str = ", ".join(map(str, matchup))
        print(f"Matchup {i + 1}: ({player_list_str})")

    print("\nIn total, the largest number of matchups that can be created is an equation which results in:")
    print(f"{len(final_matchups)}")

solve_tennis_matchups()
<<<8>>>