import math
from collections import deque

def solve_e8_torsion_count():
    """
    This script calculates the number of torsion elements of order 10 and minimal length
    in the quotient group A/Z for the Artin group A of type E8.

    The method is based on established results from Artin group theory, which reduce
    the problem to a combinatorial search on the E8 Dynkin diagram.
    """

    # The E8 Dynkin diagram can be represented as a graph.
    # Nodes are numbered 1 to 8.
    #       2
    #       |
    # 1--3--4--5--6--7--8
    e8_adj = {
        1: {3},
        2: {4},
        3: {1, 4},
        4: {2, 3, 5},
        5: {4, 6},
        6: {5, 7},
        7: {6, 8},
        8: {7}
    }
    all_nodes = set(range(1, 9))

    # Step 1: Find all subgraphs of type A4 (paths of length 3, i.e., 4 nodes).
    # We can do this with a Breadth-First Search (BFS) starting from each node.
    a4_subgraphs = set()

    for start_node in all_nodes:
        # A queue for BFS: storing (current_path_list)
        q = deque([[start_node]])
        while q:
            path = q.popleft()
            if len(path) == 4:
                # Found a path of length 3 (4 nodes).
                # Add the frozenset of nodes to avoid duplicates from different orderings.
                a4_subgraphs.add(frozenset(path))
                continue

            last_node = path[-1]
            for neighbor in e8_adj[last_node]:
                if neighbor not in path:
                    new_path = list(path)
                    new_path.append(neighbor)
                    q.append(new_path)

    # Step 2: For each A4 subgraph, find the number of commuting generators.
    # These form A4 x A1 subsystems.
    num_a4_x_a1_systems = 0
    for i, s in enumerate(a4_subgraphs):
        a4_nodes = set(s)
        
        # Find all neighbors of the A4 subgraph nodes
        neighbors_of_subgraph = set()
        for node in a4_nodes:
            neighbors_of_subgraph.update(e8_adj[node])
        
        # The commuting nodes are those not in the subgraph and not its neighbors.
        commuting_nodes = all_nodes - a4_nodes - neighbors_of_subgraph
        num_a4_x_a1_systems += len(commuting_nodes)

    # Step 3: For each A4 x A1 subsystem, there are 4! distinct Coxeter elements in the A4 part.
    # The number of elements of the required form is (# of systems) * (# of Coxeter elements in A4).
    num_coxeter_elements_in_a4 = math.factorial(4)

    total_elements = num_a4_x_a1_systems * num_coxeter_elements_in_a4
    
    print("The problem is to find the number of torsion elements of order 10 with minimal positive word length in A(E8)/Z.")
    print("This reduces to a combinatorial problem on the E8 Dynkin diagram.")
    print("The minimal length is 5, and such elements correspond to A4 x A1 subsystems.")
    print(f"1. Number of A4 subgraphs found: {len(a4_subgraphs)}")
    print(f"2. Total number of A4 x A1 subsystems: {num_a4_x_a1_systems}")
    print(f"3. Number of distinct Coxeter elements in A(A4): 4! = {num_coxeter_elements_in_a4}")
    print("\nFinal Calculation:")
    print(f"Total number of elements = (Number of A4 x A1 subsystems) * (Number of A4 Coxeter elements)")
    print(f"Total number of elements = {num_a4_x_a1_systems} * {num_coxeter_elements_in_a4} = {total_elements}")


solve_e8_torsion_count()