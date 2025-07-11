import itertools
import collections
# The cinemagoer library is required. If not installed, run: pip install cinemagoer
try:
    from cinemagoer import Cinemagoer
except ImportError:
    print("The 'cinemagoer' library is not installed.")
    print("Please install it using: pip install cinemagoer")
    exit()

def parse_years(year_string):
    """Parses a year string like '2019' or '2019-2021' into start and end years."""
    try:
        if '-' in year_string:
            start, end = map(int, year_string.split('-'))
            return start, end
        else:
            year = int(year_string)
            return year, year
    except (ValueError, TypeError):
        return None, None

def check_year_overlap(start_year, end_year, target_start=2017, target_end=2022):
    """Checks if the show's year range overlaps with the target range [2017, 2022]."""
    if start_year is None or end_year is None:
        return False
    # The condition for overlap is:
    # (ShowStart <= TargetEnd) and (ShowEnd >= TargetStart)
    return start_year <= target_end and end_year >= target_start

def analyze_graph(graph):
    """Analyzes the graph for connectivity and cyclicity."""
    nodes = list(graph.keys())
    if not nodes:
        return "Graph is empty", "Graph is empty"

    # --- Connectivity Check (using BFS) ---
    q = collections.deque([nodes[0]])
    visited_conn = {nodes[0]}
    while q:
        node = q.popleft()
        for neighbor in graph[node]:
            if neighbor not in visited_conn:
                visited_conn.add(neighbor)
                q.append(neighbor)
    
    is_connected = len(visited_conn) == len(nodes)
    connectivity = "Connected" if is_connected else "Disconnected"

    # --- Cyclicity Check (using DFS) ---
    visited_cycle = set()
    has_cycle = False
    for node in nodes:
        if node not in visited_cycle:
            # The recursion stack for the current DFS traversal
            recursion_stack = {node}
            # Use a stack for iterative DFS to avoid recursion depth limits
            stack = [(node, None)] # (current_node, parent_node)
            
            path_visited = {node}

            while stack:
                curr, parent = stack.pop()
                visited_cycle.add(curr)
                
                is_cycle_found = False
                for neighbor in graph[curr]:
                    if neighbor == parent:
                        continue
                    if neighbor in path_visited:
                        is_cycle_found = True
                        break
                    path_visited.add(neighbor)
                    stack.append((neighbor, curr))
                
                if is_cycle_found:
                    has_cycle = True
                    break
        if has_cycle:
            break
            
    cyclicity = "Cyclic" if has_cycle else "Acyclic"

    return connectivity, cyclicity

def main():
    """Main function to build and analyze the graph."""
    actors = [
        "Aaron Ashmore", "Krysten Ritter", "Emilia Jones",
        "Charlie Cox", "Devery Jacobs", "Thomas Elms"
    ]
    
    print("Initializing Cinemagoer...")
    try:
        ia = Cinemagoer()
    except Exception as e:
        print(f"Failed to connect to the data source: {e}")
        print("Please check your internet connection.")
        return

    # Cache for person and movie objects to reduce network requests
    person_cache = {}
    movie_cache = {}
    
    print("Fetching data for all actors...")
    for name in actors:
        try:
            person_results = ia.search_person(name)
            if person_results:
                person_cache[name] = ia.get_person(person_results[0].personID)
                print(f" -> Found data for {name}")
            else:
                print(f" -> Could not find data for {name}")
        except Exception as e:
            print(f"An error occurred while fetching data for {name}: {e}")
            return

    graph = collections.defaultdict(set)
    edges = []

    print("\nChecking all pairs for collaborations from 2017-2022...")
    for name1, name2 in itertools.combinations(actors, 2):
        if name1 not in person_cache or name2 not in person_cache:
            continue

        p1 = person_cache[name1]
        p2 = person_cache[name2]

        filmography1 = p1.get('filmography', {}).get('actor', [])
        filmography2 = p2.get('filmography', {}).get('actor', [])

        if not filmography1 or not filmography2:
            continue
            
        # Get sets of movie IDs for faster intersection
        ids1 = {m.movieID for m in filmography1}
        ids2 = {m.movieID for m in filmography2}
        common_ids = ids1.intersection(ids2)

        for mid in common_ids:
            try:
                if mid in movie_cache:
                    movie = movie_cache[mid]
                else:
                    movie = ia.get_movie(mid)
                    movie_cache[mid] = movie
                
                kind = movie.get('kind')
                if kind in ['tv series', 'tv mini series']:
                    year_str = movie.get('series years') or str(movie.get('year', ''))
                    start_year, end_year = parse_years(year_str)
                    
                    if check_year_overlap(start_year, end_year):
                        edges.append(tuple(sorted((name1, name2))))
                        graph[name1].add(name2)
                        graph[name2].add(name1)
                        print(f"  - Found edge: ({name1}, {name2}) via '{movie.get('title')}' ({year_str})")
                        # Once one connection is found, no need to check other common works
                        break
            except Exception as e:
                # Some movies might fail to fetch, just skip them
                # print(f"Warning: Could not fetch details for movie ID {mid}. Error: {e}")
                continue
    
    # Ensure all actors are keys in the graph, even if isolated
    for name in actors:
        if name not in graph:
            graph[name] = set()

    print("\n--- Graph Structure ---")
    if not edges:
        print("No edges were found. The graph is composed of isolated nodes.")
    else:
        print("The following connections (edges) were found:")
        # Use a set to print unique edges
        for edge in sorted(list(set(edges))):
            print(f"{edge[0]} -- {edge[1]}")
    
    print("\n--- Graph Analysis ---")
    connectivity, cyclicity = analyze_graph(graph)
    print(f"Connectivity: The graph is {connectivity}.")
    print(f"Cyclicity: The graph is {cyclicity}.")

    print("\n--- Conclusion ---")
    if connectivity == "Disconnected" and cyclicity == "Acyclic":
        print("The graph is Disconnected and acyclic, which corresponds to Choice A.")
    elif connectivity == "Disconnected" and cyclicity == "Cyclic":
        print("The graph is Disconnected and cyclic, which corresponds to Choice B.")
    elif connectivity == "Connected" and cyclicity == "Acyclic":
        print("The graph is Connected and acyclic, which corresponds to Choice C.")
    elif connectivity == "Connected" and cyclicity == "Cyclic":
        print("The graph is Connected and cyclic, which corresponds to Choice D or E.")
    else:
        print("Could not determine the graph type.")

if __name__ == '__main__':
    main()
<<<A>>>