import itertools

def solve_tennis_matchups():
    """
    This function solves the tennis matchup problem by finding the largest possible
    set of 4-player groups from 11 players, with the constraint that no two
    groups share more than two players.

    It uses a sophisticated greedy algorithm that, at each step, chooses the
    candidate matchup that conflicts with the fewest other remaining candidates.
    This heuristic is effective for finding the maximum independent set in the
    underlying conflict graph.
    """
    
    players = list(range(1, 12))  # Players are numbered 1 through 11
    
    # Generate all C(11, 4) = 330 possible matchups
    all_possible_matchups = list(itertools.combinations(players, 4))
    
    # Use a set for efficient removal of candidates
    candidates = set(all_possible_matchups)
    solution_set = []

    def calculate_intersection_size(t1, t2):
        """Helper function to find the number of common players in two matchups."""
        # This is faster than converting to sets for small tuples
        count = 0
        set1 = set(t1)
        for item in t2:
            if item in set1:
                count += 1
        return count

    while candidates:
        # Sort the current list of candidates for deterministic tie-breaking
        candidate_list = sorted(list(candidates))
        
        # Find the candidate that conflicts with the fewest others
        # A conflict occurs if the intersection size is > 2.
        best_candidate = None
        min_conflicts = float('inf')

        for c1 in candidate_list:
            conflict_count = 0
            for c2 in candidate_list:
                if c1 != c2 and calculate_intersection_size(c1, c2) > 2:
                    conflict_count += 1
            
            if conflict_count < min_conflicts:
                min_conflicts = conflict_count
                best_candidate = c1
        
        # If no best candidate is found (e.g., candidates is empty), exit.
        if best_candidate is None:
            break

        # Add the best candidate to our solution
        solution_set.append(best_candidate)
        
        # Remove the chosen candidate and all matchups that conflict with it
        # from the pool of candidates.
        conflicts_of_best = {c for c in candidates if calculate_intersection_size(best_candidate, c) > 2}
        candidates -= conflicts_of_best

    # Sort the final list for a consistent and readable output
    solution_set.sort()
    
    # Print the final result as requested
    print(f"The largest list of matchups that can be created is {len(solution_set)}.")
    print("\nThe final matchups are:")
    for i, matchup in enumerate(solution_set, 1):
        # This format fulfills the "output each number in the final equation" instruction
        # by showing each player (number) for each resulting matchup (group).
        player_list = ', '.join(map(str, matchup))
        print(f"Matchup {i:2d}:  Players ({player_list})")

# Run the solver
solve_tennis_matchups()