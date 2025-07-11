import collections

def solve_disjoint_cycles():
    """
    This function solves the DisjointCycles problem for a sample graph.
    It demonstrates a backtracking approach, which is not an FPT algorithm but is
    conceptually simpler and works for small inputs.
    """
    # Example Input
    # This graph contains two disjoint C4 cycles and one C3 cycle.
    # C4_1: 0-1-2-3-0
    # C4_2: 4-5-6-7-4
    # C3_1: 8-9-10-8
    adj = {
        0: [1, 3], 1: [0, 2], 2: [1, 3], 3: [0, 2],
        4: [5, 7], 5: [4, 6], 6: [5, 7], 7: [4, 6],
        8: [9, 10], 9: [8, 10], 10: [8, 9]
    }
    k = 2 # Parameter k

    print(f"Problem: Find k={k} vertex-disjoint simple cycles of length at least k={k}.")

    def find_all_simple_cycles(graph, min_length):
        """Finds all simple cycles of at least a given length."""
        cycles = []
        # A modified DFS to find cycles
        for start_node in graph:
            path = [start_node]
            visited = {start_node}
            
            q = collections.deque([(start_node, path, visited)])
            
            while q:
                curr_node, curr_path, curr_visited = q.popleft()
                
                for neighbor in graph.get(curr_node, []):
                    if neighbor == start_node:
                        if len(curr_path) >= min_length:
                            # To avoid duplicates, sort the cycle and store as a tuple
                            sorted_cycle = tuple(sorted(curr_path))
                            if sorted_cycle not in cycles_set:
                                cycles.append(curr_path)
                                cycles_set.add(sorted_cycle)
                    elif neighbor not in curr_visited:
                        new_path = curr_path + [neighbor]
                        new_visited = curr_visited | {neighbor}
                        q.append((neighbor, new_path, new_visited))

        cycles_set = set()
        all_cycles = []
        for start_node in adj:
            q = collections.deque([(start_node, [start_node])])
            while q:
                curr, path = q.pop()
                visited_in_path = set(path)
                for neighbor in adj.get(curr, []):
                    if neighbor == start_node:
                        if len(path) >= min_length:
                            sorted_cycle = tuple(sorted(path))
                            if sorted_cycle not in cycles_set:
                                all_cycles.append(path)
                                cycles_set.add(sorted_cycle)
                    elif neighbor not in visited_in_path:
                        q.append((neighbor, path + [neighbor]))
        return all_cycles

    # 1. Find all simple cycles of length >= k
    candidate_cycles = find_all_simple_cycles(adj, k)
    if not candidate_cycles:
        print(0)
        return

    # 2. Backtracking search for k disjoint cycles
    def backtrack(start_index, count, used_vertices):
        if count == k:
            return True
        if start_index >= len(candidate_cycles):
            return False

        for i in range(start_index, len(candidate_cycles)):
            cycle = candidate_cycles[i]
            cycle_vertices = set(cycle)
            
            # Check for disjointness
            if not used_vertices.intersection(cycle_vertices):
                # Recurse
                if backtrack(i + 1, count + 1, used_vertices.union(cycle_vertices)):
                    return True
        
        return False

    if backtrack(0, 0, set()):
        print(1)
    else:
        print(0)

solve_disjoint_cycles()