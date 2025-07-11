import collections

def solve_graph_problem():
    """
    This function builds and analyzes a graph based on actor collaborations
    in TV series between 2017 and 2022.
    """
    # Step 1: Define actors and their relevant filmography within the timeframe.
    # Data is based on filmography databases like IMDb.
    # Key: Actor, Value: Set of TV series/seasons they acted in from 2017-2022.
    filmography = {
        "Aaron Ashmore": {"Killjoys (2017-2019)", "Locke & Key (2020-2022)"},
        "Krysten Ritter": {"The Defenders (2017)", "Jessica Jones (2018-2019)"},
        "Emilia Jones": {"Locke & Key (2020-2022)"},
        "Charlie Cox": {"The Defenders (2017)", "Daredevil (2018)", "Kin (2021)", "Treason (2022)"},
        "Devery Jacobs": {"The Order (2019-2020)", "Reservation Dogs (2021-2022)"},
        "Thomas Elms": {"The Order (2019-2020)"}
    }
    nodes = list(filmography.keys())

    # Step 2: Build the graph's adjacency list and find edges
    adj_list = collections.defaultdict(list)
    edges = []

    print("Building graph for the following six nodes:")
    print(", ".join(nodes))
    print("\nAn edge exists if two actors shared a TV series/season released between 2017-2022.")
    print("\n--- Found Edges ---")

    # Iterate through all unique pairs of actors
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            actor1 = nodes[i]
            actor2 = nodes[j]

            # Find common TV shows
            common_shows = filmography[actor1].intersection(filmography[actor2])

            if common_shows:
                # Add an edge between the two actors
                adj_list[actor1].append(actor2)
                adj_list[actor2].append(actor1)
                edges.append(tuple(sorted((actor1, actor2))))
                print(f"Edge: ({actor1} - {actor2}), Reason: Shared series '{next(iter(common_shows))}'")

    if not edges:
        print("No edges found.")

    # Step 3: Analyze the graph properties

    # Connectivity check
    is_connected = True
    if not nodes:
        is_connected = False
    else:
        q = collections.deque([nodes[0]])
        visited = {nodes[0]}
        while q:
            node = q.popleft()
            for neighbor in adj_list.get(node, []):
                if neighbor not in visited:
                    visited.add(neighbor)
                    q.append(neighbor)
        if len(visited) != len(nodes):
            is_connected = False

    # Cycle check
    def has_cycle_util(u, visited_cycle, parent, current_adj_list):
        visited_cycle.add(u)
        for v in current_adj_list.get(u, []):
            if v not in visited_cycle:
                if has_cycle_util(v, visited_cycle, u, current_adj_list):
                    return True
            elif v != parent:
                return True
        return False

    is_cyclic = False
    visited_nodes_for_cycle_check = set()
    for node in nodes:
        if node not in visited_nodes_for_cycle_check:
            if has_cycle_util(node, visited_nodes_for_cycle_check, None, adj_list):
                is_cyclic = True
                break

    # Step 4: Output the conclusion
    print("\n--- Graph Analysis ---")
    
    if is_connected:
        connectivity_status = "Connected"
    else:
        connectivity_status = "Disconnected"
    
    if is_cyclic:
        cyclicity_status = "Cyclic"
    else:
        cyclicity_status = "Acyclic"

    print(f"Connectivity: The graph is {connectivity_status}.")
    print(f"Cyclicity: The graph is {cyclicity_status}.")

    print("\n--- Final Conclusion ---")
    if not is_connected and not is_cyclic:
        print("The graph is Disconnected and acyclic.")
        final_answer = "A"
    elif not is_connected and is_cyclic:
        print("The graph is Disconnected and cyclic.")
        final_answer = "B"
    elif is_connected and not is_cyclic:
        print("The graph is Connected and acyclic (a tree).")
        final_answer = "C"
    else: # Connected and cyclic
        is_cycle_graph = len(edges) == len(nodes) and all(len(adj_list.get(n, [])) == 2 for n in nodes)
        if is_cycle_graph:
            print("The graph is a cycle graph.")
            final_answer = "E"
        else:
            print("The graph is Connected and cyclic, but not a cycle graph.")
            final_answer = "D"
            
    print(f"This corresponds to Answer Choice: {final_answer}")


solve_graph_problem()