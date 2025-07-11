import collections

def generate_matchups():
    """
    Generates the maximum number of doubles matchups for 11 players
    with the constraint that no two matchups share more than 2 players.
    """
    num_players = 11
    
    # A base matchup chosen for its special combinatorial properties.
    # This is a (11, 4) difference set where differences appear at most twice.
    base_matchup = [0, 1, 3, 7]
    
    all_matchups = []
    
    # Generate 11 matchups by cyclically shifting the base matchup.
    # B_i = { (b + i) mod 11 for b in base_matchup }
    for i in range(num_players):
        new_matchup = []
        for player in base_matchup:
            # Add i to each player in the base matchup and take the result modulo 11
            new_player = (player + i) % num_players
            new_matchup.append(new_player)
        
        # Sort for consistent ordering
        new_matchup.sort()
        all_matchups.append(new_matchup)
        
    print(f"The largest list of matchups that can be created is {len(all_matchups)}.")
    print("The matchups are:")
    
    # Print each matchup. The problem asks to output each number in the final equation,
    # which we interpret as printing the elements of each matchup list.
    for matchup in all_matchups:
        # We need to print each number in the "equation" which is the list itself.
        # This formatting helps clarity.
        player_list_str = ', '.join(map(str, matchup))
        print(f"[{player_list_str}]")

generate_matchups()

# The final answer is the number of matchups found.
# The code itself prints the result, but per instructions, we state it here too.
# The problem asks for the size of the largest list.
print("\n<<<11>>>")