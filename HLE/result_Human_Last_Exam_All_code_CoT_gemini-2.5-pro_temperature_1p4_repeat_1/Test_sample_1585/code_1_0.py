import sys

def solve():
    """
    Solves the graph theory problem by determining n₂ and n₃.
    """

    # Step 1 & 2: Analyze the 2-Planar Graph and determine n₂.
    # The problem defines a 2-planar graph G₂ as having its vertices on two orthogonal planes,
    # specified as P₁ (the xy-plane, z=0) and P₂ (the yz-plane, x=0).
    # A crucial addendum states that all edges must be straight line segments that lie
    # entirely in the region where x > 0 and z > 0.

    # Let's consider an edge between two vertices, u and v, both lying on the same plane, for instance P₁ (z=0).
    # The straight line segment connecting u and v would also lie in the plane z=0.
    # This contradicts the condition that the edge must be in the region z > 0.
    # Similarly, an edge between two vertices in P₂ (x=0) would have x=0, contradicting the x > 0 requirement.
    # Therefore, edges can only exist between a vertex on P₁ and a vertex on P₂.
    # This means that G₂ must be a bipartite graph.

    # A fundamental property of bipartite graphs is that they do not contain any odd-length cycles.
    # However, the problem states that G₂ must contain exactly n induced cycles of length 5 (C₅).
    # A C₅ is an odd-length cycle.
    # This presents a direct contradiction for any graph G₂ where n > 0.
    # The only way to resolve this contradiction is if the graph has no C₅ cycles, which implies n = 0.

    # We must check if n=0 (the empty graph) is a valid solution.
    # - It is vacuously 4-regular.
    # - It contains n=0 induced C₅s.
    # - The condition on vertex removal is vacuously true.
    # - The minimality condition is vacuously true as there are no edges.
    # - The 2-planarity conditions are also vacuously satisfied.
    # Thus, the only possible value for n is 0. The smallest possible value is therefore 0.
    n2 = 0

    # Step 3 & 4: Analyze the 3-Planar Graph and determine n₃.
    # For the 3-planar graph G₃, vertices are on three mutually orthogonal planes.
    # Similar reasoning shows that the graph must be tripartite (no edges between vertices on the same plane).
    # A tripartite graph can have odd cycles, including a C₅. So, a non-zero n₃ might exist.
    # However, the question asks for the *smallest possible value* of n.
    # Just as with the 2-planar case, we can check if n=0 is a valid solution for G₃.
    # The empty graph (n=0) vacuously satisfies all the listed properties for 3-planarity as well.
    # Since n=0 is a possible solution, and no smaller non-negative integer exists, it must be the smallest value.
    n3 = 0

    # Step 5: Calculate the final result.
    result = (n2 + n3) * n2

    # Print the reasoning and the final calculation.
    print("Step 1: Finding n_2 for the 2-planar graph.")
    print("The geometric constraints for a 2-planar graph require it to be bipartite.")
    print("A bipartite graph cannot contain odd cycles, such as a C5.")
    print("The problem requires n C5 cycles. This is only possible if n = 0.")
    print("The empty graph (n=0) vacuously satisfies all the conditions.")
    print(f"Therefore, the smallest possible value is n_2 = {n2}")
    print("-" * 20)
    
    print("Step 2: Finding n_3 for the 3-planar graph.")
    print("The problem asks for the smallest possible n. The empty graph (n=0) also vacuously satisfies all conditions for the 3-planar case.")
    print("Since 0 is a possible value, it must be the smallest.")
    print(f"Therefore, the smallest possible value is n_3 = {n3}")
    print("-" * 20)

    print("Step 3: Calculating the final result.")
    print(f"The required calculation is (n_2 + n_3) * n_2")
    print(f"Plugging in the values: ({n2} + {n3}) * {n2} = {result}")

solve()
