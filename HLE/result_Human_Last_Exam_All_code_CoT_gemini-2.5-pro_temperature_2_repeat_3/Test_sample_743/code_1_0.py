import itertools

def find_max_matchups():
    """
    Finds the largest possible list of 4-player matchups from 11 players,
    such that no two matchups have more than two players in common.
    """
    num_players = 11
    group_size = 4
    players = range(num_players)

    # 1. Generate all possible matchups
    all_possible_matchups = list(itertools.combinations(players, group_size))

    # 2. Use a greedy approach to build the solution list
    solution_list = []
    for candidate_matchup in all_possible_matchups:
        candidate_set = set(candidate_matchup)
        is_valid = True

        # 3. Check against all matchups already in our solution
        for existing_matchup in solution_list:
            existing_set = set(existing_matchup)
            
            # 4. The constraint: intersection size must not be more than 2
            if len(candidate_set.intersection(existing_set)) > 2:
                is_valid = False
                break
        
        # 5. If valid, add the candidate to our solution
        if is_valid:
            solution_list.append(candidate_matchup)

    # 6. Print the final result
    print(f"The largest list of matchups that can be created is: {len(solution_list)}")
    print("The list of matchups is:")
    for i, matchup in enumerate(solution_list):
        # Convert tuple of player numbers to a list of strings for display
        player_list = [f"Player {p}" for p in matchup]
        print(f"Matchup {i+1}: {list(matchup)}")

if __name__ == '__main__':
    find_max_matchups()
