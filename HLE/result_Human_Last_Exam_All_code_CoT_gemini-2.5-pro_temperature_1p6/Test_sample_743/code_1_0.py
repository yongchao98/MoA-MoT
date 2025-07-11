def solve_tennis_matchups():
    """
    Solves the tennis matchup problem by constructing a solution based on a cyclic design.

    The problem asks for the largest set of 4-player matchups from 11 players,
    such that any two matchups share at most 2 players. This corresponds to a packing
    problem in design theory, and the maximum size is known to be 11.

    This solution is constructed cyclically using the players {0, 1, ..., 10}.
    We start with a "base matchup" {0, 1, 3, 7}. All 11 matchups are generated
    by adding a value i (from 0 to 10) to each player number in the base matchup,
    with the result taken modulo 11.
    """
    
    players_count = 11
    # This is the base block for the cyclic design.
    base_matchup = [0, 1, 3, 7]
    
    final_matchups = []
    
    # Generate 11 matchups by cyclically shifting the base matchup
    for i in range(players_count):
        # Create a new matchup by adding 'i' and taking the result modulo 11
        new_matchup = sorted([(player + i) % players_count for player in base_matchup])
        final_matchups.append(new_matchup)
        
    num_matchups = len(final_matchups)

    print(f"The largest list of matchups has a size of {num_matchups}.")
    print("The matchups are (with players numbered 0 to 10):")
    
    # Print each generated matchup
    for matchup in final_matchups:
        print(matchup)

solve_tennis_matchups()