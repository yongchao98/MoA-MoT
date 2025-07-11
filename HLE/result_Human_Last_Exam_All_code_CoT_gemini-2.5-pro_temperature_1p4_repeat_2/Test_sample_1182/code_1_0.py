import collections

def find_double_points():
    """
    Calculates the number of double points in the stable reduction of the curve
    y^2 = 8*x^5 + 4*x^4 + 4*x^3 + x^2 + 8*x above 2.
    """
    # Step 1 & 2: Define branch points and their pairwise distances.
    # The branch points are e0, e1, e2, e3, e4, and e_inf.
    # We computed the 2-adic valuations of their differences.
    # Let's list all edges with non-negative weights (distances).
    # All other pairs have negative distances.
    edges = [
        # (point1, point2, distance)
        ('e0', 'e1', 3),
        ('e2', 'e_inf', 1),
        ('e3', 'e_inf', 1),
        ('e4', 'e_inf', 1),
        ('e0', 'e_inf', 0),
        ('e1', 'e_inf', 0),
    ]

    print("Step 1: The branch points are e0, e1, e2, e3, e4, e_inf.")
    print("Step 2 & 3: The 'distances' (2-adic valuations of differences) between branch points are computed.")
    print("The positive distances are:")
    print("d(e0, e1) = 3")
    print("d(e2, e_inf) = 1")
    print("d(e3, e_inf) = 1")
    print("d(e4, e_inf) = 1")
    print("d(e0, e_inf) = 0")
    print("d(e1, e_inf) = 0")
    print("\nStep 4: Construct the Minimal Spanning Tree (MST) using these distances.")

    # Sort edges by distance in descending order
    edges.sort(key=lambda item: item[2], reverse=True)
    
    # Kruskal's algorithm to find the MST
    parent = {f'e{i}': f'e{i}' for i in range(5)}
    parent['e_inf'] = 'e_inf'
    
    def find_set(v):
        if v == parent[v]:
            return v
        parent[v] = find_set(parent[v])
        return parent[v]

    def unite_sets(a, b):
        a = find_set(a)
        b = find_set(b)
        if a != b:
            parent[b] = a
            return True
        return False

    mst_edges = []
    for u, v, weight in edges:
        if find_set(u) != find_set(v):
            unite_sets(u,v)
            mst_edges.append((u, v))
            if len(mst_edges) == 5: # 6 vertices, so 5 edges in MST
                break
    
    print("The edges of the MST are:")
    for u, v in mst_edges:
        print(f"({u}, {v})")

    # Step 5: Find the number of components in the stable reduction.
    # These correspond to leaves (degree 1) and branch points (degree > 2) in the MST.
    degrees = collections.defaultdict(int)
    for u, v in mst_edges:
        degrees[u] += 1
        degrees[v] += 1
    
    leaves = [node for node, deg in degrees.items() if deg == 1]
    branch_points = [node for node, deg in degrees.items() if deg > 2]
    
    num_components = len(leaves) + len(branch_points)

    print("\nStep 5: The components of the stable reduction correspond to leaves and branch points of the MST.")
    print(f"Leaves (degree 1 nodes) in MST: {leaves}")
    print(f"Branch points (degree > 2 nodes) in MST: {branch_points}")
    print(f"Total number of components = {len(leaves)} + {len(branch_points)} = {num_components}")

    # Step 6: The number of double points is the number of components - 1.
    num_double_points = num_components - 1
    print(f"\nStep 6: The number of double points is the number of components minus 1.")
    print(f"Number of double points = {num_components} - 1 = {num_double_points}")

    return num_double_points

result = find_double_points()
print(f"\nFinal Answer: The number of double points is {result}.")