import collections

# Before running, please install the necessary library:
# pip install cinemagoer

from cinemagoer import Cinemagoer

def solve_graph_problem():
    """
    This function constructs and analyzes the actor collaboration graph.
    """
    actors = {
        'Aaron Ashmore': '0039149',
        'Krysten Ritter': '1269983',
        'Emilia Jones': '4317869',
        'Charlie Cox': '1214452',
        'Devery Jacobs': '2711235',
        'Thomas Elms': '7261972'
    }
    actor_names = list(actors.keys())
    
    # Initialize Cinemagoer instance to fetch IMDb data
    ia = Cinemagoer()
    
    # Store actor filmographies to avoid redundant API calls
    filmographies = {}
    print("Fetching filmography data for each actor...")
    for name, person_id in actors.items():
        person = ia.get_person(person_id)
        # Focus on acting roles
        if 'actor' in person.data['filmography']:
             filmographies[name] = person.data['filmography']['actor']
        elif 'actress' in person.data['filmography']:
             filmographies[name] = person.data['filmography']['actress']
        else:
             filmographies[name] = []
        print(f"- Fetched data for {name}")

    # Adjacency list to represent the graph
    adj = collections.defaultdict(set)
    edges = set()

    print("\nFinding edges by checking shared TV series (2017-2022):")
    
    # Iterate through all unique pairs of actors
    for i in range(len(actor_names)):
        for j in range(i + 1, len(actor_names)):
            actor1_name = actor_names[i]
            actor2_name = actor_names[j]

            # Create a set of valid TV series for the first actor for quick lookups
            actor1_works = set()
            for work in filmographies[actor1_name]:
                kind = work.get('kind', 'N/A')
                year = work.get('year')
                if (kind in ('tv series', 'tv mini series')) and year and (2017 <= year <= 2022):
                    actor1_works.add(work.movieID)

            # Check if the second actor has any work in common
            for work in filmographies[actor2_name]:
                if work.movieID in actor1_works:
                    # Found a valid collaboration
                    edge = tuple(sorted((actor1_name, actor2_name)))
                    if edge not in edges:
                        edges.add(edge)
                        adj[actor1_name].add(actor2_name)
                        adj[actor2_name].add(actor1_name)
                        print(f"- Edge found: {actor1_name} -- {actor2_name} (via '{work['title']}' [{work['year']}])")

    if not edges:
        print("- No edges found in the graph.")

    # --- Graph Analysis ---
    print("\nAnalyzing the graph...")

    # 1. Connectivity Check
    visited = set()
    if not actor_names:
        is_connected = True
    else:
        q = collections.deque([actor_names[0]])
        visited.add(actor_names[0])
        while q:
            node = q.popleft()
            for neighbor in adj[node]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    q.append(neighbor)
        
        is_connected = (len(visited) == len(actor_names))

    if is_connected:
        print("Result 1: The graph is connected.")
    else:
        print(f"Result 1: The graph is disconnected. (Visited {len(visited)} of {len(actor_names)} nodes from one component)")

    # 2. Cycle Check using DFS
    visited_for_cycle = set()
    has_cycle = False
    
    def detect_cycle_util(u, parent):
        nonlocal has_cycle
        if has_cycle: return

        visited_for_cycle.add(u)
        
        for v in adj[u]:
            if v == parent:
                continue
            if v in visited_for_cycle:
                has_cycle = True
                return
            detect_cycle_util(v, u)

    for node in actor_names:
        if node not in visited_for_cycle:
            detect_cycle_util(node, None)
            if has_cycle:
                break
    
    if has_cycle:
        print("Result 2: The graph is cyclic.")
    else:
        print("Result 2: The graph is acyclic.")

    # --- Final Conclusion ---
    print("\nConclusion:")
    if is_connected:
        if has_cycle:
            # Check if it's a cycle graph (every node has degree 2)
            is_cycle_graph = all(len(adj[node]) == 2 for node in actor_names)
            if is_cycle_graph:
                print("The graph is a cycle graph.")
                final_answer = "E"
            else:
                print("The graph is connected and cyclic, but not a cycle graph.")
                final_answer = "D"
        else:
            print("The graph is connected and acyclic (a tree).")
            final_answer = "C"
    else: # Disconnected
        if has_cycle:
            print("The graph is disconnected and cyclic.")
            final_answer = "B"
        else:
            print("The graph is disconnected and acyclic.")
            final_answer = "A"

    print(f"\nThe correct description is: {final_answer}")
    return final_answer

if __name__ == '__main__':
    final_answer_choice = solve_graph_problem()
    # In a real scenario, this would be the final output requested by the platform
    # print(f"<<<{final_answer_choice}>>>")
<<<A>>>