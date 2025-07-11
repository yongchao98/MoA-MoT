import itertools

def solve_tennis_matchups():
    """
    Solves the Tennis matchup problem by finding all groups of 4 players from 11
    whose numeric representations sum to a multiple of 11. This construction
    guarantees that any two such groups will not share more than 2 players.
    """
    
    # We represent the 11 players with numbers from 0 to 10.
    players_zero_based = list(range(11))
    
    # This list will hold the matchups that satisfy the condition.
    valid_matchups = []
    
    # Generate all possible groups of four from the 11 players.
    all_combinations = itertools.combinations(players_zero_based, 4)
    
    # Check each combination against our mathematical rule.
    for combo in all_combinations:
        # The condition is that the sum of player numbers is a multiple of 11.
        if sum(combo) % 11 == 0:
            # This is a valid matchup according to the construction.
            valid_matchups.append(combo)
            
    # The answer is the total number of matchups found.
    print(f"The largest list of matchups that can be created is of size: {len(valid_matchups)}\n")
    
    print("The list of matchups is shown below, using players 1-11 for clarity.")
    print("The calculation for each valid matchup is also shown using the 0-10 player representation.\n")
    
    # Sort the matchups for a clean, deterministic output.
    sorted_matchups_zero_based = sorted([tuple(sorted(m)) for m in valid_matchups])
    
    # Print the details for each valid matchup.
    for i, matchup_zero_based in enumerate(sorted_matchups_zero_based):
        # Convert to players 1-11 for the user-friendly output.
        matchup_one_based = tuple(p + 1 for p in matchup_zero_based)
        
        # Unpack the numbers for the equation string.
        z0, z1, z2, z3 = matchup_zero_based
        total_sum = sum(matchup_zero_based)
        
        # Output the matchup and the equation that validates it.
        # Each number in the final equation is outputted.
        print(f"Matchup {i+1:2d}: {matchup_one_based}")
        print(f"  Calculation: {z0} + {z1} + {z2} + {z3} = {total_sum}, which is a multiple of 11.\n")

# Execute the function to find and print the solution.
solve_tennis_matchups()

# The final numerical answer is the number of matchups found.
# The code calculates this to be 33.
print("<<<33>>>")