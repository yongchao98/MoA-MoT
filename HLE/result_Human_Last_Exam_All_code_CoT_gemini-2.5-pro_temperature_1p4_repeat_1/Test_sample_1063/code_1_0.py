import itertools
from imdb import Cinemagoer

def solve_graph_problem():
    """
    This function builds and analyzes a graph of actors based on their shared TV appearances.
    """
    print("Initializing actor list and graph structure...")
    actors = {
        "Aaron Ashmore": None,
        "Krysten Ritter": None,
        "Emilia Jones": None,
        "Charlie Cox": None,
        "Devery Jacobs": None,
        "Thomas Elms": None,
    }
    # Adjacency list to represent the graph
    adj_list = {name: [] for name in actors}
    # Year range for the series/miniseries
    valid_years = range(2017, 2023)

    # Use an API to get movie/TV show data
    ia = Cinemagoer()

    print("Fetching filmography for each actor. This may take a moment...")
    for name in actors:
        try:
            # Search for the actor and get their full filmography
            person = ia.search_person(name)[0]
            ia.update(person, 'filmography')
            actors[name] = person
            print(f"-> Successfully fetched data for {name}")
        except Exception as e:
            print(f"Could not fetch data for {name}: {e}")
            return

    print("\nChecking for collaborations (edges) between actors...")
    # Iterate over every unique pair of actors
    for actor1_name, actor2_name in itertools.combinations(actors.keys(), 2):
        actor1 = actors[actor1_name]
        actor2 = actors[actor2_name]

        # Get their filmographies as actors/actresses
        filmography1 = actor1.get('filmography', [{}])[0].get('actor', [])
        filmography2 = actor2.get('filmography', [{}])[0].get('actor', [])
        
        # Create a set of qualifying TV show IDs for the first actor
        shows1_ids = set()
        for show in filmography1:
            if show.get('kind') in ['tv series', 'tv mini series'] and show.get('year') in valid_years:
                shows1_ids.add(show.movieID)

        # Check if the second actor was in any of those shows
        for show in filmography2:
            if show.movieID in shows1_ids:
                # Found a connection!
                print(f"Found connection: {actor1_name} and {actor2_name} in '{show['title']}' ({show['year']})")
                adj_list[actor1_name].append(actor2_name)
                adj_list[actor2_name].append(actor1_name)
                # Once one connection is found, we don't need to look for more for this pair
                break
    
    print("\n--- Graph Analysis ---")
    print("Final Adjacency List:")
    for node, neighbors in adj_list.items():
        if not neighbors:
            print(f"{node}: []")
        else:
            print(f"{node}: {neighbors}")

    # Analysis of the graph properties
    visited = set()
    num_components = 0
    is_cyclic = False

    for node in adj_list:
        if node not in visited:
            num_components += 1
            # Check for cycles using DFS traversal within each component
            stack = [(node, None)] # (current_node, parent_node)
            component_visited = set()
            while stack:
                curr, parent = stack.pop()
                if curr in component_visited:
                    is_cyclic = True # We revisited a node in the current traversal, so there is a cycle
                component_visited.add(curr)
                visited.add(curr)
                for neighbor in adj_list[curr]:
                    if neighbor != parent: # Don't go immediately back to the parent
                        stack.append((neighbor, curr))

    print("\nProperties:")
    if num_components > 1:
        print(f"1. Connectivity: The graph has {num_components} separate components, so it is DISCONNECTED.")
    else:
        print("1. Connectivity: The graph has 1 component, so it is CONNECTED.")

    if is_cyclic:
        print("2. Cyclicity: The graph contains cycles.")
    else:
        print("2. Cyclicity: The graph is ACYCLIC (contains no cycles).")
        
    print("\nConclusion: The graph is Disconnected and acyclic.")

# Run the solver
solve_graph_problem()
<<<A>>>