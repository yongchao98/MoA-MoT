import itertools

def find_max_matchups():
    """
    This function calculates the largest list of tennis matchups according to the specified rules.
    It uses a greedy algorithm to find a maximal set of compatible matchups.
    """
    
    # There are 11 players, represented by numbers 0 to 10.
    players = list(range(11))
    
    # A matchup is a group of 4 players. We generate all possible combinations.
    all_possible_matchups = list(itertools.combinations(players, 4))
    
    # This list will store our final selection of matchups.
    final_matchups = []
    
    # We iterate through every possible matchup and try to add it to our final list.
    for candidate_matchup in all_possible_matchups:
        # Assume the candidate is compatible until proven otherwise.
        is_compatible = True
        
        # Check the candidate against every matchup already in our final list.
        for existing_matchup in final_matchups:
            # The restriction: no two matchups share more than two players.
            # We find the number of common players using set intersection.
            common_players = set(candidate_matchup).intersection(set(existing_matchup))
            
            if len(common_players) > 2:
                # The candidate shares too many players with an existing matchup.
                # It's not compatible, so we reject it and stop checking.
                is_compatible = False
                break
        
        # If, after checking all existing matchups, the candidate remained compatible...
        if is_compatible:
            # ...we add it to our final list.
            final_matchups.append(candidate_matchup)
            
    # After checking all possibilities, we have our result.
    print(f"The largest list of matchups that can be created is: {len(final_matchups)}")
    print("The matchups in this list are:")
    
    # Print each matchup in the final list.
    for matchup in final_matchups:
        # The numbers represent the players in each matchup.
        print(list(matchup))

# Run the function to get the answer.
find_max_matchups()