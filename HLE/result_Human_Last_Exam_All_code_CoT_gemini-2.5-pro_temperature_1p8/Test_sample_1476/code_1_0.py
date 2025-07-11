import numpy as np

def explain_and_demonstrate():
    """
    Explains the logical deduction and demonstrates it with a sample graph.
    """
    print("Thinking Process and Explanation:")
    print("=================================")
    print("1. From the problem statement, we are given three key conditions:")
    print("   a) 'no cycles having non-zero sum': This implies the edge signal x^1 is a gradient of a vertex potential z. Mathematically: x^1 = B₁ᵀ z.")
    print("   b) 'B₁ x¹ 1ᵀ = 0': This simplifies to B₁ x¹ = 0. This means x^1 is divergence-free.")
    print("   c) 'x¹_e = |x⁰_u - x⁰_v|': The edge signal is the absolute difference of the vertex signal x⁰.")
    print("\n2. Combining (a) and (b), we get B₁(B₁ᵀ z) = 0.")
    print("   The term B₁B₁ᵀ is the graph Laplacian L₀. So, L₀ z = 0.")
    print("   This implies that the potential z must be constant on each connected component of the graph.")
    print("\n3. If z is constant on a component, its gradient B₁ᵀ z must be zero for all edges within that component.")
    print("   Therefore, from (a), the edge signal x¹ must be the zero vector: x¹ = 0.")
    print("\n4. Using (c), since x¹ = 0, we have |x⁰_u - x⁰_v| = 0 for every edge {u,v}.")
    print("   This means the difference in the original vertex signal x⁰ is zero across every edge.")
    print("\n5. The Total Variation (TV) of the graph is defined as TV = Σ |x⁰_u - x⁰_v| over all edges.")
    print("   Since each term is 0, the total variation must be 0.")
    print("   This corresponds to option D.")
    
    print("\n\nNumerical Demonstration:")
    print("========================")
    # Let's use a simple path graph with 4 nodes and 3 edges as an example.
    # V = {0, 1, 2, 3}, E = {{0,1}, {1,2}, {2,3}}
    # Edges oriented as 0->1, 1->2, 2->3
    
    # Vertex-edge incidence matrix B1 (|V|x|E|)
    B1 = np.array([[-1., 0., 0.],
                   [ 1.,-1., 0.],
                   [ 0., 1.,-1.],
                   [ 0., 0., 1.]])
    
    print("For a path graph with 4 nodes, the B₁ matrix is:\n", B1)
    
    # The graph Laplacian L0 = B1 * B1^T
    L0 = B1 @ B1.T
    print("\nThe corresponding graph Laplacian L₀ = B₁B₁ᵀ is:\n", L0)
    
    print("\nThe condition L₀ z = 0 implies z is constant (e.g., z = [c, c, c, c]ᵀ).")
    c = 5 # Arbitrary constant
    z = np.array([c, c, c, c])
    print(f"Let's test with z = {z}, L₀z = {L0 @ z}")
    
    print("\nNow we find x¹ from x¹ = B₁ᵀ z:")
    x1 = B1.T @ z
    print(f"x¹ = B₁ᵀ z = {x1}")
    
    print("\nAs shown, the premises force the edge signal x¹ to be the zero vector.")
    
    # From x¹_e = |x⁰_u - x⁰_v|, we deduce the components for the total variation calculation.
    # For our path graph, these components are |x⁰_0 - x⁰_1|, |x⁰_1 - x⁰_2|, |x⁰_2 - x⁰_3|.
    # Since x¹ = [0, 0, 0], we have:
    diff_e1 = 0  # |x⁰_0 - x⁰_1|
    diff_e2 = 0  # |x⁰_1 - x⁰_2|
    diff_e3 = 0  # |x⁰_2 - x⁰_3|

    print("\nFinally, we calculate the Total Variation:")
    tv = diff_e1 + diff_e2 + diff_e3
    
    # Print the equation as requested
    print(f"TV = |x⁰₀-x⁰₁| + |x⁰₁-x⁰₂| + |x⁰₂-x⁰₃| = {diff_e1} + {diff_e2} + {diff_e3} = {tv}")
    print("\nConclusion: The total variation of the graph is 0.")

explain_and_demonstrate()