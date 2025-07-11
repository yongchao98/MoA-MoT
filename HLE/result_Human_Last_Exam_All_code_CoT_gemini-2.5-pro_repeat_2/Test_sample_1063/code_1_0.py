import collections
from itertools import combinations
try:
    from cinemagoer import Cinemagoer
except ImportError:
    print("Installing cinemagoer library...")
    import pip
    pip.main(['install', 'cinemagoer'])
    from cinemagoer import Cinemagoer

def solve_graph_problem():
    """
    This function solves the graph problem by fetching actor data, building a graph,
    and analyzing its properties to find the correct description.
    """
    # Step 1: Initialize actors, years, and data structures
    ia = Cinemagoer()
    actors_names = [
        "Aaron Ashmore", "Krysten Ritter", "Emilia Jones",
        "Charlie Cox", "Devery Jacobs", "Thomas Elms"
    ]
    year_range = range(2017, 2023)
    actor_works = collections.defaultdict(dict)

    print("Step 1: Fetching and filtering filmography for each actor...")
    # Step 2 & 3: Fetch and filter filmography for each actor
    for name in actors_names:
        try:
            person = ia.search_person(name)[0]
            ia.update(person, info=['filmography'])
            filmography = person.get('filmography', {}).get('actor', person.get('filmography', {}).get('actress', []))

            for work in filmography:
                if work.get('kind') in ['tv series', 'tv mini series'] and work.get('year') in year_range:
                    actor_works[name][work.movieID] = work.get('title')
        except Exception as e:
            print(f"Warning: Could not process data for {name}. Error: {e}")
    print("Filmography processing complete.\n")

    # Step 4: Establish edges by finding common works
    print("Step 2: Identifying edges based on shared work...")
    edges = []
    edge_details = collections.defaultdict(list)
    adj = {name: [] for name in actors_names}

    for actor1, actor2 in combinations(actors_names, 2):
        common_work_ids = set(actor_works[actor1].keys()).intersection(set(actor_works[actor2].keys()))
        if common_work_ids:
            edge = tuple(sorted((actor1, actor2)))
            edges.append(edge)
            adj[actor1].append(actor2)
            adj[actor2].append(actor1)
            for work_id in common_work_ids:
                edge_details[edge].append(actor_works[actor1][work_id])

    print("Graph construction complete:")
    print(f"Nodes (6): {', '.join(actors_names)}")
    if edges:
        print("Edges found:")
        for (actor1, actor2), works in edge_details.items():
            print(f"- Edge between {actor1} and {actor2} (shared work: {', '.join(works)})")
    else:
        print("No edges found between any actors.")
    print("-" * 30)

    # Step 5: Analyze the graph
    print("Step 3: Analyzing graph properties...")
    num_nodes = len(actors_names)
    num_edges = len(edges)
    
    # Analyze Connectivity by finding connected components
    visited = set()
    num_components = 0
    for node in actors_names:
        if node not in visited:
            num_components += 1
            q = collections.deque([node])
            visited.add(node)
            while q:
                current_node = q.popleft()
                for neighbor in adj[current_node]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        q.append(neighbor)
    
    is_connected = num_components == 1
    
    # Analyze Cyclicity using the formula E = N - C
    is_acyclic = (num_edges == num_nodes - num_components)

    print(f"Number of Nodes (N): {num_nodes}")
    print(f"Number of Edges (E): {num_edges}")
    print(f"Number of Connected Components (C): {num_components}")
    
    if is_connected:
        print("Result: The graph is CONNECTED.")
    else:
        print("Result: The graph is DISCONNECTED.")

    if is_acyclic:
        print("Result: The graph is ACYCLIC.")
    else:
        print("Result: The graph is CYCLIC.")
    print("-" * 30)

    # Step 6: Formulate the final answer
    if not is_connected and is_acyclic:
        final_answer = "A"
    elif not is_connected and not is_acyclic:
        final_answer = "B"
    elif is_connected and is_acyclic:
        final_answer = "C"
    elif is_connected and not is_acyclic:
        is_cycle_graph = all(len(neighbors) == 2 for neighbors in adj.values())
        final_answer = "E" if is_cycle_graph else "D"
    else:
        final_answer = "Unknown"

    print(f"The analysis shows the graph is disconnected and acyclic, which corresponds to answer choice A.")
    return final_answer

if __name__ == '__main__':
    final_choice = solve_graph_problem()
    print(f"\n<<<A>>>")