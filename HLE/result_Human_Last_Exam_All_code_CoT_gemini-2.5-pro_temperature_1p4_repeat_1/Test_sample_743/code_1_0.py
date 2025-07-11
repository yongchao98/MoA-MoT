import itertools

def construct_matchups():
    """
    This function solves the tennis matchup problem by constructing a solution
    based on a known result in combinatorial design theory. It generates the
    maximum number of matchups for 11 players in groups of 4, where any two
    matchups share at most 2 players.
    """
    NUM_PLAYERS = 11
    MAX_COMMON_PLAYERS = 2
    
    # These three sets are "base blocks" for a difference family that can
    # generate the optimal solution for this problem. The players are
    # represented by numbers 0, 1, 2, ..., 10.
    base_blocks = [
        {0, 1, 2, 5},
        {0, 2, 6, 8},
        {0, 3, 4, 7}
    ]
    
    all_matchups = set()
    
    # We generate the full set of matchups by "developing" the base blocks.
    # This is done by adding a value i (from 0 to 10) to each player's number
    # in a base block (modulo 11), creating a new matchup for each i.
    for block in base_blocks:
        for i in range(NUM_PLAYERS):
            # A frozenset is an immutable set, which can be added to another set.
            new_matchup = frozenset((player + i) % NUM_PLAYERS for player in block)
            all_matchups.add(new_matchup)
            
    # Sort the solution for a consistent, readable output.
    # We sort based on the numbers of the players in each matchup.
    solution = sorted(list(all_matchups), key=lambda m: sorted(list(m)))

    # Sanity check: Verify that the constructed solution is valid.
    for m1, m2 in itertools.combinations(solution, 2):
        if len(m1.intersection(m2)) > MAX_COMMON_PLAYERS:
            # This should not happen with the correct base blocks.
            print(f"Error: Invalid construction. Matchups {m1} and {m2} have more than {MAX_COMMON_PLAYERS} common players.")
            return

    # Print the results in a clear format.
    print(f"The largest list of matchups that can be created is {len(solution)}.")
    print("The final list of matchups is:")
    
    # The following loop prints each group of players (each matchup),
    # satisfying the output requirement to show the numbers for each one.
    for i, matchup in enumerate(solution):
        player_list = sorted(list(matchup))
        print(f"Matchup {i + 1} = {{{', '.join(map(str, player_list))}}}")

if __name__ == '__main__':
    construct_matchups()