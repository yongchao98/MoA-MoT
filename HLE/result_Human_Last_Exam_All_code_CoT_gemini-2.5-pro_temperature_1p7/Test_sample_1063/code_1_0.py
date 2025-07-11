import collections
import itertools

def solve_graph_problem():
    """
    This function solves the graph problem by:
    1. Defining the actors and their relevant filmographies.
    2. Building the graph by finding common TV shows from 2017-2022.
    3. Printing the edges found and the reasons for them.
    4. Analyzing the graph's properties (connectivity and cycles).
    5. Determining the correct answer choice.
    """
    actors = [
        "Aaron Ashmore", "Krysten Ritter", "Emilia Jones",
        "Charlie Cox", "Devery Jacobs", "Thomas Elms"
    ]
    num_nodes = len(actors)

    # Filmography data: Actor -> set of relevant (Show, Year) tuples
    # Years are for the first episode of the miniseries or season.
    filmography = {
        "Aaron Ashmore": {("Locke & Key", 2020)},
        "Krysten Ritter": {("The Defenders", 2017)},
        "Emilia Jones": {("Locke & Key", 2020)},
        "Charlie Cox": {("The Defenders", 2017)},
        "Devery Jacobs": {("The Order", 2019)},
        "Thomas Elms": {("The Order", 2019)},
    }

    print("Step 1: Identifying edges based on shared TV series (2017-2022)\n")

    adjacency_list = collections.defaultdict(list)
    edges = []
    
    # Iterate through all unique pairs of actors
    for actor1, actor2 in itertools.combinations(actors, 2):
        # Find common shows using set intersection on the show tuples
        common_shows = filmography[actor1].intersection(filmography[actor2])
        if common_shows:
            show, year = list(common_shows)[0]
            edge = tuple(sorted((actor1, actor2)))
            edges.append(edge)
            adjacency_list[actor1].append(actor2)
            adjacency_list[actor2].append(actor1)
            print(f"Edge Found: ({actor1}, {actor2})")
            print(f"Reason: Both acted in '{show}' which began in {year}.\n")

    num_edges = len(edges)
    
    print("\nStep 2: Analyzing graph properties")
    print(f"Number of nodes: {num_nodes}")
    print(f"Number of edges: {num_edges}\n")

    # Analyze connectivity by finding connected components using Breadth-First Search (BFS)
    visited = set()
    num_components = 0
    for node in actors:
        if node not in visited:
            num_components += 1
            q = collections.deque([node])
            visited.add(node)
            while q:
                current_node = q.popleft()
                for neighbor in adjacency_list[current_node]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        q.append(neighbor)
    
    is_connected = num_components == 1
    
    # Analyze cycles. A graph is acyclic if #edges = #nodes - #components
    is_acyclic = (num_edges == num_nodes - num_components)
    
    print("Analysis Result:")
    print(f"The graph has {num_components} connected component(s).")
    print(f"Is the graph connected? {'Yes' if is_connected else 'No'}")
    print(f"Is the graph acyclic? {'Yes' if is_acyclic else 'No'}\n")

    # Determine the final answer
    if not is_connected and is_acyclic:
        result = "A. Disconnected and acyclic"
        final_answer = "<<<A>>>"
    elif not is_connected and not is_acyclic:
        result = "B. Disconnected and cyclic"
        final_answer = "<<<B>>>"
    elif is_connected and is_acyclic:
        result = "C. Connected and acyclic"
        final_answer = "<<<C>>>"
    elif is_connected and not is_acyclic:
        # Check if it is a cycle graph (all vertices degree 2)
        is_cycle_graph = all(len(adj) == 2 for adj in adjacency_list.values())
        if is_cycle_graph:
            result = "E. A cycle graph"
            final_answer = "<<<E>>>"
        else:
            result = "D. Connected and cyclic, but not a cycle graph"
            final_answer = "<<<D>>>"

    print("Conclusion:")
    print(f"The graph is correctly described as: {result}")
    print(final_answer)

solve_graph_problem()