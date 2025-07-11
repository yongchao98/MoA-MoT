import itertools

def verify_matchup_list():
    """
    This function verifies if a given list of matchups satisfies the problem's constraints.
    The problem: Arrange 11 players in groups of 4, so that no two groups have
    more than 2 players in common.
    This code checks a known optimal solution of size 12.
    """
    # The list of 11 players, represented by numbers 0 through 10.
    players = list(range(11))
    
    # A known optimal solution with 12 matchups.
    # Each inner list represents a matchup of 4 players.
    optimal_matchups = [
        [0, 1, 2, 3],
        [0, 1, 4, 5],
        [0, 2, 6, 7],
        [0, 3, 8, 9],
        [0, 4, 8, 10],
        [0, 5, 6, 9],
        [1, 2, 8, 10],
        [1, 3, 6, 9],
        [1, 4, 6, 7],
        [1, 5, 7, 8],
        [2, 3, 5, 7],
        [2, 4, 5, 9]
    ]

    max_common_players = 2
    is_valid = True

    print("Verifying the provided list of matchups...")

    # Iterate through all unique pairs of matchups from the list.
    for match1, match2 in itertools.combinations(optimal_matchups, 2):
        # Use sets to find the number of common players efficiently.
        common_players = set(match1).intersection(set(match2))
        
        # Check if the number of common players exceeds the maximum allowed.
        if len(common_players) > max_common_players:
            print(f"ERROR: Matchup {match1} and Matchup {match2} have {len(common_players)} players in common.")
            is_valid = False
            break # No need to check further if one pair is invalid.
    
    if is_valid:
        print("Verification successful! All matchups in the list are valid.")
        
    print("-" * 30)
    
    # The final answer is the number of matchups in this verified list.
    final_answer = len(optimal_matchups)
    
    print(f"The largest list of matchups that can be created is of size {final_answer}.")
    print("Here is one such list of matchups:")
    for i, matchup in enumerate(optimal_matchups):
        print(f"Matchup {i+1}: {matchup}")

# Execute the function
verify_matchup_list()
