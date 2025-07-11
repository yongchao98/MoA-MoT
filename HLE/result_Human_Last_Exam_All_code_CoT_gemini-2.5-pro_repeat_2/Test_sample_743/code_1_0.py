import itertools

def generate_matchups():
    """
    This function generates the largest possible set of doubles matchups
    from 11 players, such that no two matchups share more than 2 players.
    
    The method uses a cyclic construction based on modular arithmetic.
    """
    num_players = 11
    
    # We define a "base" matchup using players numbered 0-10 for calculation ease.
    # This specific set {0, 1, 2, 4} is chosen because it generates a valid design
    # when cyclically shifted.
    base_matchup = {0, 1, 2, 4}
    
    matchups = []
    for i in range(num_players):
        # Generate a new matchup by adding 'i' to each player in the base matchup,
        # and taking the result modulo num_players. This 'rotates' the players.
        new_matchup = sorted([(p + i) % num_players for p in base_matchup])
        matchups.append(tuple(new_matchup))
        
    # Remove duplicate matchups if any (in this construction, there won't be any).
    # Convert to 1-based indexing for user-friendly output.
    final_matchups = sorted([tuple(p + 1 for p in m) for m in set(matchups)])
    
    print(f"The largest list of matchups that can be created is {len(final_matchups)}.")
    print("Here is one such list of matchups:")
    
    # Print each matchup in the final "equation" or list
    for i, matchup in enumerate(final_matchups):
        print(f"Matchup {i+1}: {matchup}")

# Run the function to print the result.
generate_matchups()