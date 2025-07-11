import collections

def is_simple_cycle(graph, path):
    """Checks if a given path is a simple cycle in the graph."""
    if len(path) < 3:
        return False
    
    # Check for simplicity (no repeated vertices)
    if len(set(path)) != len(path):
        return False

    # Check for connectivity
    for i in range(len(path) - 1):
        if path[i+1] not in graph.get(path[i], []):
            return False
    
    # Check if it's a closed loop
    if path[0] not in graph.get(path[-1], []):
        return False
        
    return True

def check_disjoint_cycles(graph, k, cycles):
    """
    Verifies if the given 'cycles' list is a valid solution for the DisjointCycles problem.
    A valid solution consists of k vertex-disjoint simple cycles, each of length at least k.
    """
    # 1. Check if there are k cycles
    if len(cycles) != k:
        print(f"Failed: Expected {k} cycles, but found {len(cycles)}.")
        return False

    all_vertices = set()
    total_vertices_count = 0
    
    for i, cycle_path in enumerate(cycles):
        # 2. Check if each is a valid simple cycle
        if not is_simple_cycle(graph, cycle_path):
            print(f"Failed: Candidate {i} ({cycle_path}) is not a simple cycle.")
            return False
            
        # 3. Check if cycle length is at least k
        if len(cycle_path) < k:
            print(f"Failed: Cycle {i} ({cycle_path}) has length {len(cycle_path)}, which is less than k={k}.")
            return False
            
        # Add vertices to check for disjointness later
        for vertex in cycle_path:
            all_vertices.add(vertex)
        total_vertices_count += len(cycle_path)

    # 4. Check if cycles are vertex-disjoint
    if len(all_vertices) != total_vertices_count:
        print("Failed: Cycles are not vertex-disjoint.")
        return False

    print("Success: The provided cycles are a valid solution.")
    return True

if __name__ == '__main__':
    # Parameter for the problem
    k = 3

    # A sample graph represented as an adjacency list.
    # This graph contains 3 disjoint cycles of length >= 3.
    # C1: 0-1-2-0 (length 3)
    # C2: 3-4-5-3 (length 3)
    # C3: 6-7-8-9-6 (length 4)
    graph = {
        0: [1, 2], 1: [0, 2], 2: [0, 1],
        3: [4, 5], 4: [3, 5], 5: [3, 4],
        6: [7, 9], 7: [6, 8], 8: [7, 9], 9: [6, 8],
        10: [0] # extra node/edge
    }

    print(f"Problem: Find k={k} vertex-disjoint simple cycles of length at least k={k}.")
    print(f"Graph has vertices: {list(graph.keys())}\n")

    # --- Test Cases ---

    # Case 1: A valid solution
    valid_solution = [[0, 1, 2], [3, 4, 5], [6, 7, 8, 9]]
    print(f"Checking candidate 1: {valid_solution}")
    check_disjoint_cycles(graph, k, valid_solution)
    print("-" * 20)

    # Case 2: An invalid solution (not enough cycles)
    invalid_solution_1 = [[0, 1, 2], [3, 4, 5]]
    print(f"Checking candidate 2: {invalid_solution_1}")
    check_disjoint_cycles(graph, k, invalid_solution_1)
    print("-" * 20)

    # Case 3: An invalid solution (cycles not disjoint)
    # Let's create a graph where cycles might overlap
    graph_2 = {
        0: [1, 5], 1: [0, 2], 2: [1, 3], 3: [2, 4], 4: [3, 5], 5: [0, 4, 6],
        6: [5, 7], 7: [6, 8], 8: [7, 5]
    }
    # C1: 0-1-2-3-4-5-0 (len 6)
    # C2: 5-6-7-8-5 (len 4)
    # These share vertex 5
    invalid_solution_2 = [[0, 1, 2, 3, 4, 5], [5, 6, 7, 8]]
    print(f"Checking candidate 3: {invalid_solution_2} with k=2 in a different graph")
    check_disjoint_cycles(graph_2, 2, invalid_solution_2)
    print("-" * 20)
    
    # Case 4: An invalid solution (one path is not a cycle)
    invalid_solution_3 = [[0, 1, 2], [3, 4, 5], [6, 7]]
    print(f"Checking candidate 4: {invalid_solution_3}")
    check_disjoint_cycles(graph, k, invalid_solution_3)
    print("-" * 20)
    
    # Case 5: An invalid solution (cycle length too short)
    k_short = 4
    invalid_solution_4 = [[0, 1, 2], [3, 4, 5], [6, 7, 8, 9]]
    print(f"Checking candidate 5: {invalid_solution_4} with k={k_short}")
    check_disjoint_cycles(graph, k_short, invalid_solution_4)
    print("-" * 20)
