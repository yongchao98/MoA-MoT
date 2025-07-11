def explain_graph_problem():
    """
    Explains the reasoning behind the properties of a maximal planar graph
    when a new edge is added.
    """
    print("Analyzing the graph G' = G U {e}, where G is a maximal planar graph.")
    print("="*70)

    print("Step 1: Properties of a Maximal Planar Graph (G)")
    print("- A planar graph can be drawn on a plane with no edges crossing.")
    print("- 'Maximal planar' means no more edges can be added without making it non-planar.")
    print("- In any planar drawing of a maximal planar graph, all faces are triangles.")
    print("-" * 70)

    print("Step 2: Adding an edge 'e' to G")
    print("- We form G' = G U {e}, where e = (u, v) is an edge not already in G.")
    print("- Because G is maximal planar, G' must be non-planar.")
    print("- Any curve connecting u and v must cross at least one existing edge of G.")
    print("- Therefore, the crossing number of G' is at least 1.")
    print("- This rules out options claiming G' has a plane drawing (no crossings).")
    print("-" * 70)

    print("Step 3: The Crossing Number of G'")
    print("- Maximal planar graphs on n>=4 vertices are 3-connected.")
    print("- A theorem states: If H is a 3-connected planar graph, adding one new edge gives a graph with a crossing number of exactly 1.")
    print("- Applying this theorem to G, the crossing number of G' is 1.")
    print("- This means G' can be drawn with at most one crossing, ruling out the option that it might require more.")
    print("-" * 70)

    print("Step 4: Uniqueness of the 1-Crossing Drawing")
    print("- The final question is whether this 1-crossing drawing is unique.")
    print("- A drawing's topology is defined by which pair of edges cross.")
    print("- Consider G = K5 - e' (K5 with one edge removed), which is maximal planar.")
    print("- Let's add that edge e' back, so G' = K5, which has a crossing number of 1.")
    print("- We can draw G' with one crossing in several ways. For example, by placing one vertex inside a triangle formed by three others, and the fifth vertex outside.")
    print("- The edge connecting the inside and outside vertices must cross the triangle.")
    print("- We can choose which of the three edges of the triangle it crosses. Each choice creates a different drawing.")
    print("- Since there are multiple ways to form a 1-crossing drawing, the drawing is not unique.")
    print("="*70)

    print("Conclusion:")
    print("G' can be drawn with at most one crossing, and this drawing is not unique.")

if __name__ == '__main__':
    explain_graph_problem()