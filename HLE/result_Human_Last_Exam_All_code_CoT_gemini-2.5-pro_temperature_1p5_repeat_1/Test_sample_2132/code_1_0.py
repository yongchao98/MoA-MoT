import collections

def calculate_bottcher_complexity():
    """
    Calculates the Böttcher Molecular Complexity for cyclopentanecarboxylic acid.
    """
    # Step 1: Define the molecular graph for cyclopentanecarboxylic acid.
    # We represent the 8 non-hydrogen atoms with numbers:
    # 1-5: Carbons in the cyclopentane ring (C1-C5)
    # 6:   The carboxyl carbon (C6), attached to C1
    # 7-8: The two oxygen atoms of the carboxyl group, attached to C6
    graph = collections.defaultdict(list, {
        1: [2, 5, 6],
        2: [1, 3],
        3: [2, 4],
        4: [3, 5],
        5: [4, 1],
        6: [1, 7, 8],
        7: [6],
        8: [6]
    })

    # Step 2: Calculate N_atoms and N_bonds.
    nodes = sorted(list(graph.keys()))
    n_atoms = len(nodes)
    
    # N_bonds is the number of edges, which is half the sum of degrees
    n_bonds = sum(len(adj) for adj in graph.values()) // 2
    
    print(f"Product of the reaction: Cyclopentanecarboxylic acid")
    print(f"Number of non-hydrogen atoms (N_atoms): {n_atoms}")
    print(f"Number of bonds between non-hydrogen atoms (N_bonds): {n_bonds}")

    # Step 3: Calculate N_paths.
    # This function finds all simple paths between a start and end node.
    def find_all_paths(current_graph, start, end, path=[]):
        path = path + [start]
        if start == end:
            return [path]
        if start not in current_graph:
            return []
        paths = []
        for node in current_graph[start]:
            if node not in path: # Prevents cycles
                new_paths = find_all_paths(current_graph, node, end, path)
                for p in new_paths:
                    paths.append(p)
        return paths

    n_paths = 0
    # Iterate through all unique pairs of nodes to sum up the path counts.
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            start_node = nodes[i]
            end_node = nodes[j]
            paths_between_pair = find_all_paths(graph, start_node, end_node)
            n_paths += len(paths_between_pair)
            
    print(f"Total number of simple paths (N_paths): {n_paths}")
    
    # Step 4: Calculate the Böttcher Molecular Complexity.
    if n_paths == 0:
        complexity = float('inf')
    else:
        complexity = (n_atoms * n_bonds) / n_paths
    
    print("\n--- Calculation ---")
    print(f"Böttcher Molecular Complexity = (N_atoms * N_bonds) / N_paths")
    print(f"Result = ({n_atoms} * {n_bonds}) / {n_paths}")
    print(f"Final Value = {complexity}")

if __name__ == "__main__":
    calculate_bottcher_complexity()