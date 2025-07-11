import itertools

def solve_tennis_matchups():
    """
    Finds the largest set of 4-player matchups from 11 players such that
    no two matchups share more than two players.
    """
    num_players = 11
    group_size = 4
    max_common_players = 2

    players = range(num_players)
    
    # Generate all possible 4-player groups from the 11 players
    all_possible_matchups = list(itertools.combinations(players, group_size))
    
    # This will hold the final list of matchups that satisfy the condition
    selected_matchups = []

    # Implement the greedy algorithm
    for candidate_matchup in all_possible_matchups:
        candidate_set = set(candidate_matchup)
        is_valid = True
        
        # Check the candidate against all previously selected matchups
        for existing_matchup in selected_matchups:
            existing_set = set(existing_matchup)
            
            # If the intersection is too large, the candidate is invalid
            if len(candidate_set.intersection(existing_set)) > max_common_players:
                is_valid = False
                break
        
        # If the candidate is valid, add it to our list
        if is_valid:
            selected_matchups.append(candidate_matchup)

    # Print the results
    print("Found the following list of matchups:")
    for i, matchup in enumerate(selected_matchups):
        # The final equation is the matchup itself, showing the players.
        # Example: 0 + 1 + 2 + 3 is one matchup.
        equation_str = " + ".join(map(str, matchup))
        print(f"Matchup {i+1:2d}: Players {equation_str}")
        
    print("\n-----------------------------------------------------")
    print(f"The largest list of matchups that can be created is: {len(selected_matchups)}")

# Execute the function to solve the problem
solve_tennis_matchups()

<<<11>>>