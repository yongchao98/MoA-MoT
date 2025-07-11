import itertools

def find_largest_matchup_list():
    """
    Finds the largest list of 4-player matchups from 11 players,
    ensuring no two matchups have more than two players in common.
    """
    players = list(range(1, 12))
    group_size = 4
    # The constraint |A intersect B| <= 2 means no two matchups can share a 3-player group.
    # So, we track all 3-player combinations.
    common_players_group_size = 3

    # 1. Generate all possible 4-player matchups
    all_possible_matchups = list(itertools.combinations(players, group_size))

    # 2. These will store our results
    final_matchups = []
    # Using a set for efficient lookup of used 3-player groups (triples)
    used_triples = set()

    # 3. Greedily iterate through all potential matchups
    for matchup in all_possible_matchups:
        # Generate all 3-player triples for the current matchup
        triples_in_matchup = list(itertools.combinations(matchup, common_players_group_size))

        # 4. Check if this matchup is valid
        is_valid = True
        for triple in triples_in_matchup:
            if triple in used_triples:
                is_valid = False
                break
        
        # 5. If it's valid, add it to our list and update used triples
        if is_valid:
            final_matchups.append(matchup)
            used_triples.update(triples_in_matchup)

    # 6. Print the final results
    print(f"The largest list of matchups found contains {len(final_matchups)} matchups.")
    print("-----------------------------------------")
    print("The matchups are:")
    for m in final_matchups:
        # The prompt requires printing each number in the "final equation"
        print(f"Matchup: {m[0]}, {m[1]}, {m[2]}, {m[3]}")

# Execute the function to find and print the solution
find_largest_matchup_list()
<<<18>>>