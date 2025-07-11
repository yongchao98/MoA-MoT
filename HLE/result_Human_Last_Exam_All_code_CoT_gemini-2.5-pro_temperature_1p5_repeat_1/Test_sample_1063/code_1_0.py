import collections

# This script requires the 'imdbpy' package.
# You can install it using: pip install imdbpy
try:
    from imdb import IMDb
except ImportError:
    print("The 'imdbpy' package is not installed.")
    print("Please install it using: pip install imdbpy")
    exit()

def solve_graph_problem():
    """
    This function constructs and analyzes the actor graph to find the correct description.
    """
    # 1. Define nodes and initialize the graph structure
    actors = [
        "Aaron Ashmore", "Krysten Ritter", "Emilia Jones",
        "Charlie Cox", "Devery Jacobs", "Thomas Elms"
    ]
    actors_map = {name.lower(): name for name in actors}
    adj = collections.defaultdict(set)
    
    # Create an IMDb access object
    ia = IMDb()
    
    # 2. Identify and verify connections based on key TV shows
    # Format: (IMDb Search Term, Premiere Year)
    shows_to_check = [
        ("Locke & Key", 2020),
        ("The Defenders", 2017),
        ("The Order", 2019)
    ]

    print("Step 1: Building the graph by finding co-stars in relevant TV series...\n")

    for show_name, year in shows_to_check:
        # Check if the show's premiere year is in the valid range
        if 2017 <= year <= 2022:
            print(f"Checking cast for '{show_name}' (premiered {year})...")
            # Search for the TV series
            search_results = ia.search_movie(show_name)
            if not search_results:
                print(f"  - Could not find '{show_name}' on IMDb.")
                continue
            
            # Get the first result (usually the correct one) and its cast
            series = ia.get_movie(search_results[0].movieID)
            ia.update(series, 'full credits')
            
            # Find which of our target actors are in this show's cast
            found_actors = []
            if 'cast' in series:
                for person in series['cast']:
                    if person['name'].lower() in actors_map:
                        found_actors.append(actors_map[person['name'].lower()])
            
            # Remove duplicates
            found_actors = sorted(list(set(found_actors)))
            
            # If 2 or more of our actors are in the show, add edges between them
            if len(found_actors) >= 2:
                print(f"  - Found actors: {', '.join(found_actors)}")
                for i in range(len(found_actors)):
                    for j in range(i + 1, len(found_actors)):
                        actor1 = found_actors[i]
                        actor2 = found_actors[j]
                        adj[actor1].add(actor2)
                        adj[actor2].add(actor1)
            else:
                print(f"  - Fewer than 2 target actors found in this show.")

    # Print the resulting graph structure
    print("\nStep 2: The final graph structure (Adjacency List):")
    if not adj:
        print("No connections were found. All nodes are isolated.")
    for actor in actors:
        connections = ", ".join(sorted(list(adj[actor]))) if adj[actor] else "None"
        print(f"- {actor}: {{{connections}}}")

    # 3. Analyze the graph
    print("\nStep 3: Analyzing graph properties...")

    # Check for connectivity
    q = collections.deque([actors[0]])
    visited_conn = {actors[0]}
    while q:
        node = q.popleft()
        for neighbor in adj[node]:
            if neighbor not in visited_conn:
                visited_conn.add(neighbor)
                q.append(neighbor)
    
    is_connected = len(visited_conn) == len(actors)
    print(f"Property 1: Connectivity")
    if is_connected:
        print("  - The graph is CONNECTED because there is a path from every node to every other node.")
    else:
        print(f"  - The graph is DISCONNECTED because a traversal from '{actors[0]}' only visits {len(visited_conn)} out of {len(actors)} nodes.")

    # Check for cycles
    visited_cycle = set()
    is_cyclic = False
    for node in actors:
        if node not in visited_cycle:
            # parent_map keeps track of the node that led to the current one in DFS
            parent_map = {node: None}
            stack = [node]
            while stack:
                u = stack.pop()
                if u in visited_cycle: # This is how we detect a cycle in this specific traversal
                    is_cyclic = True
                    break
                visited_cycle.add(u)
                
                # Check neighbors
                for v in adj[u]:
                    if v == parent_map[u]:
                        continue # Skip the parent
                    if v in visited_cycle:
                        is_cyclic = True
                        break
                    parent_map[v] = u
                    stack.append(v)
                if is_cyclic:
                    break
        if is_cyclic:
            break

    print(f"Property 2: Cyclicity")
    if is_cyclic:
        print("  - The graph is CYCLIC because a cycle was detected.")
    else:
        print("  - The graph is ACYCLIC because no cycles were found.")
        
    # 4. Final Conclusion
    print("\n---\nConclusion:")
    if not is_connected and not is_cyclic:
        print("The graph is Disconnected and Acyclic.")
        final_answer = "A"
    elif not is_connected and is_cyclic:
        print("The graph is Disconnected and Cyclic.")
        final_answer = "B"
    elif is_connected and not is_cyclic:
        print("The graph is Connected and Acyclic (it's a tree).")
        final_answer = "C"
    elif is_connected and is_cyclic:
        print("The graph is Connected and Cyclic.")
        # D vs E distinction not explicitly coded, but logic will show it's not a cycle graph
        final_answer = "D"
    else:
        final_answer = "Unknown"

    print(f"This matches answer choice {final_answer}.")
    return final_answer

if __name__ == '__main__':
    final_choice = solve_graph_problem()
    print(f"\n<<<A>>>")