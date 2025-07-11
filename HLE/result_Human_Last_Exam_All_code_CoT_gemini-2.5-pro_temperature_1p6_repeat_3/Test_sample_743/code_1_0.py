import itertools

def solve_tennis_matchups():
    """
    Finds the largest list of 4-player matchups from 11 players such that
    no two matchups share more than two players.
    """
    
    # 1. Generate all possible matchups
    players = list(range(11))
    all_matchups_list = list(itertools.combinations(players, 4))
    
    # Use tuples for hashability and a canonical ordering
    all_matchups = tuple(sorted(all_matchups_list))
    num_matchups = len(all_matchups)
    
    # Create a map from matchup to its index for quick lookups
    matchup_to_idx = {m: i for i, m in enumerate(all_matchups)}

    # 2. Build the conflict graph
    # An edge exists if two matchups share 3 players.
    conflict_adj_sets = [set() for _ in range(num_matchups)]
    for i in range(num_matchups):
        for j in range(i + 1, num_matchups):
            m1 = set(all_matchups[i])
            m2 = set(all_matchups[j])
            
            # A conflict occurs if the intersection size is 3
            if len(m1.intersection(m2)) == 3:
                conflict_adj_sets[i].add(j)
                conflict_adj_sets[j].add(i)

    # 3. Solver for Maximum Independent Set using recursion with memoization
    # The memoization table will store results for subproblems
    memo = {}

    def find_max_independent_set(candidates):
        """
        Recursively finds the maximum independent set from a set of candidate nodes.
        'candidates' is a frozenset of node indices.
        Returns the frozenset of nodes in the MIS.
        """
        # Base case: no more candidates
        if not candidates:
            return frozenset()
            
        # If we have already solved for this set of candidates, return the stored result
        if candidates in memo:
            return memo[candidates]

        # Heuristic: Pick a pivot node 'v' to branch on.
        # Choosing the one with the most connections within the candidate set is a good heuristic
        # as it helps prune the search space faster.
        pivot = -1
        max_degree = -1
        for node in candidates:
            # count neighbors of 'node' that are also in the 'candidates' set
            degree = len(conflict_adj_sets[node].intersection(candidates))
            if degree > max_degree:
                max_degree = degree
                pivot = node
        
        # All candidates are isolated, so they form an independent set
        if max_degree == 0:
            return candidates
        
        # Get neighbors of the pivot that are in the current candidate set
        pivot_neighbors = conflict_adj_sets[pivot].intersection(candidates)

        # Branch 1: Include the pivot in the independent set.
        # This means we cannot include any of its neighbors.
        remaining_candidates1 = candidates - {pivot} - pivot_neighbors
        solution_with_pivot = {pivot}.union(find_max_independent_set(remaining_candidates1))

        # Branch 2: Exclude the pivot from the independent set.
        # We then find the MIS in the remaining graph.
        remaining_candidates2 = candidates - {pivot}
        solution_without_pivot = find_max_independent_set(remaining_candidates2)

        # Choose the better of the two branches
        if len(solution_with_pivot) > len(solution_without_pivot):
            result = solution_with_pivot
        else:
            result = solution_without_pivot
            
        # Store the result in the memoization table before returning
        memo[candidates] = result
        return result

    # Initial call with all matchup indices as candidates
    all_candidate_indices = frozenset(range(num_matchups))
    max_set_indices = find_max_independent_set(all_candidate_indices)
    
    # Prepare the final result
    final_matchups = [all_matchups[i] for i in max_set_indices]
    
    # Print the results
    print("The largest list of matchups found is:")
    for i, matchup in enumerate(sorted(final_matchups)):
        # Output each number (player ID) in the matchup
        player_str = " ".join(map(str, matchup))
        print(f"Matchup {i+1:>2}: {player_str}")
        
    print("\nThis means the largest list contains the following number of matchups:")
    print(len(final_matchups))


if __name__ == "__main__":
    solve_tennis_matchups()
