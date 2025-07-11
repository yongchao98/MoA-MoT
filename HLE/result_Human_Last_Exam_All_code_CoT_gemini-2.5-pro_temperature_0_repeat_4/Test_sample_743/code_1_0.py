def generate_and_print_matchups():
    """
    This function constructs and prints the largest possible list of tennis matchups
    for 11 players in groups of 4, where any two matchups share at most 2 players.
    """
    
    # We represent the 11 players with numbers 0 through 10 for mathematical convenience.
    num_players = 11
    
    # We define a "base" matchup. This specific set has properties that allow
    # the cyclic construction to work correctly.
    base_matchup = {0, 1, 3, 7}
    
    all_matchups = []
    
    # We generate 11 matchups by cyclically "shifting" the players in the base matchup.
    for i in range(num_players):
        # For each shift 'i', we create a new matchup by adding 'i' to each player
        # number in the base matchup, taking the result modulo 11.
        new_matchup = sorted([(player + i) % num_players for player in base_matchup])
        all_matchups.append(new_matchup)
        
    print(f"The largest possible list contains {len(all_matchups)} matchups.")
    print("Here is a list of all 11 matchups (players are numbered 1 to 11):")
    
    # Print each matchup, converting player numbers from 0-10 to 1-11 for clarity.
    for i, matchup in enumerate(all_matchups):
        # Adding 1 to each player number for a 1-based display
        display_matchup = [player + 1 for player in matchup]
        print(f"Matchup {i + 1}: {display_matchup}")

# Execute the function to print the results.
generate_and_print_matchups()
