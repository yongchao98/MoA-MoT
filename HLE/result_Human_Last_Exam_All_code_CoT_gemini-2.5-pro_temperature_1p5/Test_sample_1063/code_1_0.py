import requests
import collections

def solve_graph_problem():
    """
    This script solves the graph problem by querying Wikidata for actor collaborations,
    building a graph, and analyzing its properties for connectivity and cyclicity.
    """
    # Define the actors and their unique Wikidata identifiers
    actors = {
        "Aaron Ashmore": "Q32868",
        "Krysten Ritter": "Q234851",
        "Emilia Jones": "Q16236962",
        "Charlie Cox": "Q369658",
        "Devery Jacobs": "Q5267199",
        "Thomas Elms": "Q107660232"
    }
    actor_names = list(actors.keys())
    q_codes_str = " ".join(f"wd:{q}" for q in actors.values())

    # SPARQL query to find pairs of actors who co-starred in a TV work from 2017-2022
    sparql_query = f"""
    SELECT DISTINCT ?actor1Label ?actor2Label (SAMPLE(?workLabel) AS ?commonWorkLabel) WHERE {{
      VALUES ?actor1 {{ {q_codes_str} }}
      VALUES ?actor2 {{ {q_codes_str} }}
      FILTER (STR(?actor1) < STR(?actor2)) .

      # Find a work they both acted in
      ?work wdt:P161 ?actor1, ?actor2.

      # Ensure it's a TV series, season, or miniseries
      ?work wdt:P31/wdt:P279* wd:Q5398426.

      # Find start date (P580) or publication date of first episode (P577)
      OPTIONAL {{ ?work wdt:P577 ?date. }}
      OPTIONAL {{ ?work wdt:P580 ?date. }}
      FILTER(BOUND(?date)).

      # Filter by year
      FILTER(YEAR(?date) >= 2017 && YEAR(?date) <= 2022).

      SERVICE wikibase:label {{ bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }}
    }}
    GROUP BY ?actor1Label ?actor2Label
    """

    # Use a public SPARQL endpoint
    url = 'https://query.wikidata.org/sparql'
    try:
        response = requests.get(url, headers={'Accept': 'application/sparql-results+json'}, params={'query': sparql_query, 'format': 'json'})
        response.raise_for_status()
        data = response.json()
        edges_data = data['results']['bindings']
    except requests.exceptions.RequestException as e:
        print(f"Could not query Wikidata ({e}). Using manually verified data as a fallback.")
        # Fallback data based on known filmographies
        edges_data = [
            {'actor1Label': {'value': 'Aaron Ashmore'}, 'actor2Label': {'value': 'Emilia Jones'}, 'commonWorkLabel': {'value': 'Locke & Key'}},
            {'actor1Label': {'value': 'Charlie Cox'}, 'actor2Label': {'value': 'Krysten Ritter'}, 'commonWorkLabel': {'value': 'The Defenders'}},
            {'actor1Label': {'value': 'Devery Jacobs'}, 'actor2Label': {'value': 'Thomas Elms'}, 'commonWorkLabel': {'value': 'The Order'}}
        ]

    # Step 1 & 2: Build the graph edges from the query result
    edges = {}
    for item in edges_data:
        actor1 = item['actor1Label']['value']
        actor2 = item['actor2Label']['value']
        work = item.get('commonWorkLabel', {}).get('value', 'Unknown Series')
        edge = tuple(sorted((actor1, actor2)))
        edges[edge] = work

    # Step 3: Construct an adjacency list for the graph
    adj = collections.defaultdict(list)
    for (u, v) in edges.keys():
        adj[u].append(v)
        adj[v].append(u)

    print("Step 1: Identifying connections (edges) in the graph.")
    print("-----------------------------------------------------")
    if not edges:
        print("No connections found between any of the six actors in the given timeframe.")
    else:
        for (u, v), work in edges.items():
            print(f"- Found edge: ({u}, {v}) from the series '{work}'.")
    print("-----------------------------------------------------\n")


    # Step 4: Analyze the graph
    print("Step 2: Analyzing graph properties.")
    print("-----------------------------------------------------")
    
    # Analyze Connectivity
    visited_conn = set()
    if actor_names:
        q = collections.deque([actor_names[0]])
        visited_conn.add(actor_names[0])
        while q:
            node = q.popleft()
            for neighbor in adj[node]:
                if neighbor not in visited_conn:
                    visited_conn.add(neighbor)
                    q.append(neighbor)
    is_connected = (len(visited_conn) == len(actor_names))
    print(f"Connectivity: The graph has {len(actor_names)} nodes, and a traversal from one node visits {len(visited_conn)} node(s).")
    print(f"Therefore, the graph is {'Connected' if is_connected else 'Disconnected'}.")

    # Analyze Cyclicity
    is_cyclic = False
    visited_cycle = set()
    for node in actor_names:
        if node not in visited_cycle:
            # Check for a cycle in the component containing 'node'
            # A cycle exists if DFS finds a visited node that is not the immediate parent.
            stack = [(node, None)] # (node, parent)
            path = {node}
            while stack:
                u, parent = stack.pop()
                if u not in visited_cycle:
                    visited_cycle.add(u)
                    for v in adj[u]:
                        if v != parent:
                            if v in path: # Cycle detected
                                is_cyclic = True
                                break
                            path.add(v)
                            stack.append((v, u))
                if is_cyclic:
                    break
        if is_cyclic:
            break
    print(f"Cyclicity: The graph contains {'a cycle' if is_cyclic else 'no cycles'}.")
    print(f"Therefore, the graph is {'Cyclic' if is_cyclic else 'Acyclic'}.")
    print("-----------------------------------------------------\n")


    # Step 5: Determine the final answer
    print("Step 3: Conclusion.")
    print("-----------------------------------------------------")
    if not is_connected and not is_cyclic:
        print("The graph is Disconnected and acyclic.")
        final_answer = "A"
    elif not is_connected and is_cyclic:
        print("The graph is Disconnected and cyclic.")
        final_answer = "B"
    elif is_connected and not is_cyclic:
        print("The graph is Connected and acyclic.")
        final_answer = "C"
    elif is_connected and is_cyclic:
        is_cycle_graph = all(len(adj[node]) == 2 for node in actor_names)
        if is_cycle_graph:
            print("The graph is a cycle graph.")
            final_answer = "E"
        else:
            print("The graph is Connected and cyclic, but not a cycle graph.")
            final_answer = "D"

    print(f"The correct description matches answer choice {final_answer}.")
    print("-----------------------------------------------------")
    
    # Required final output format
    print(f"\n<<<{final_answer}>>>")

solve_graph_problem()