import itertools

def solve_tennis_matchups():
    """
    Solves the tennis matchup problem by modeling it as a Maximum Independent Set problem
    on an incompatibility graph.
    """

    # Memoization cache for the recursive function
    memo = {}
    
    # Global graph to be accessible by the recursive function
    incompatibility_graph = {}

    def find_max_set_recursive(nodes_tuple):
        """
        Finds the maximum size of a compatible set of matchups from a given
        tuple of nodes (matchups) using a backtracking algorithm with memoization.
        This is equivalent to finding the maximum independent set.
        """
        if not nodes_tuple:
            return []

        # Return cached result if this subproblem was already solved
        if nodes_tuple in memo:
            return memo[nodes_tuple]

        # Select the first node to make a decision on.
        v = nodes_tuple[0]
        
        # All other nodes to consider.
        remaining_nodes_list = list(nodes_tuple[1:])
        remaining_nodes_tuple = tuple(remaining_nodes_list)

        # === Case 1: The solution does NOT include node 'v' ===
        # We find the max set from the rest of the nodes.
        solution1 = find_max_set_recursive(remaining_nodes_tuple)

        # === Case 2: The solution DOES include node 'v' ===
        # We must remove all matchups incompatible with 'v' from our consideration.
        incompatible_with_v = incompatibility_graph[v]
        
        # Create the new list of potential nodes for the recursive call
        nodes_for_case2 = tuple([node for node in remaining_nodes_list if node not in incompatible_with_v])

        # The solution for case 2 is 'v' plus the best solution from the pruned list.
        solution2_rec = find_max_set_recursive(nodes_for_case2)
        solution2 = [v] + solution2_rec

        # Compare the two cases and choose the larger set.
        if len(solution1) > len(solution2):
            result = solution1
        else:
            result = solution2
        
        # Cache the result before returning
        memo[nodes_tuple] = result
        return result

    # --- Main script execution ---

    num_players = 11
    players = list(range(num_players))

    # 1. Generate all C(11, 4) = 330 possible matchups.
    # We use frozensets because they are hashable and can be used as dictionary keys.
    all_matchups = [frozenset(m) for m in itertools.combinations(players, 4)]
    
    # 2. Build the incompatibility graph.
    nonlocal incompatibility_graph
    incompatibility_graph = {m: set() for m in all_matchups}
    
    for i in range(len(all_matchups)):
        for j in range(i + 1, len(all_matchups)):
            m1 = all_matchups[i]
            m2 = all_matchups[j]
            # Incompatible if they share 3 or more players.
            if len(m1.intersection(m2)) > 2:
                incompatibility_graph[m1].add(m2)
                incompatibility_graph[m2].add(m1)

    # 3. Find the maximum compatible set (max independent set).
    # Convert list to a sorted tuple to ensure it's hashable for memoization.
    all_matchups_tuple = tuple(sorted(list(all_matchups), key=lambda s: tuple(sorted(list(s)))))
    
    max_set = find_max_set_recursive(all_matchups_tuple)

    # 4. Print the result. The numbers in each list represent the players.
    print(f"The largest list of matchups that can be created is of size: {len(max_set)}")
    print("Here is one such list of matchups:")
    
    # Sort the final set for a clean, consistent output
    sorted_max_set = sorted([sorted(list(m)) for m in max_set])
    
    for matchup in sorted_max_set:
        print(f"  {matchup}")

if __name__ == '__main__':
    solve_tennis_matchups()