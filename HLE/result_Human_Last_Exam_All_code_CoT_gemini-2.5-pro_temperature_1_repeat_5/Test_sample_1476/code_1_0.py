import numpy as np

def demonstrate_inference():
    """
    This function demonstrates the logical deduction on a sample graph.
    We will use a triangle graph, which has a cycle, to show that option B is not necessary.
    """
    print("--- Demonstration on a Triangle Graph ---")
    
    # 1. Define the graph G=(V, E)
    # Vertices V = {0, 1, 2}
    # Edges E = {e0, e1, e2} where e0={0,1}, e1={1,2}, e2={2,0}
    # We orient the edges: 0->1, 1->2, 2->0 to define B1
    
    print("Graph: Vertices V={0, 1, 2}, Edges E={{0,1}, {1,2}, {2,0}}")
    edges = [(0, 1), (1, 2), (2, 0)]
    
    # 2. Define the vertex-edge incidence matrix B1
    # Rows are vertices (0, 1, 2), columns are edges (e0, e1, e2)
    B1 = np.array([
        [-1, 0, 1],  # Vertex 0 is tail of e0, head of e2
        [1, -1, 0],  # Vertex 1 is head of e0, tail of e1
        [0, 1, -1]   # Vertex 2 is head of e1, tail of e2
    ])
    print("\nVertex-Edge Incidence Matrix B1:\n", B1)

    # 3. Define signals x0 and x1 according to the problem statement.
    # Let's choose a simple x0 that leads to x1=0, which satisfies the premises.
    # According to our derivation, x0 must be constant on the component.
    x0 = np.array([5.0, 5.0, 5.0])
    print("\nLet's define a vertex signal x0:", x0)

    # 4. Calculate x1 based on Condition 3: x1_e = |x0_u - x0_v|
    x1 = np.array([
        abs(x0[0] - x0[1]),
        abs(x0[1] - x0[2]),
        abs(x0[2] - x0[0])
    ])
    print("Calculate edge signal x1 using x1_e = |x0_u - x0_v|:", x1)

    # 5. Verify that the premises hold for these signals.
    print("\n--- Verifying Premises ---")
    
    # Premise 1: "no cycles having non-zero sum"
    # The cycle is 0->1->2->0. The sum is x1[0] + x1[1] + x1[2].
    cycle_sum = x1[0] + x1[1] + x1[2]
    print(f"Premise 1 check: Sum of x1 over the cycle = {x1[0]} + {x1[1]} + {x1[2]} = {cycle_sum}. Premise holds.")

    # Premise 2: B1 * x1 = 0
    div_x1 = B1 @ x1
    print(f"Premise 2 check: B1 * x1 = {div_x1}. Premise holds.")

    # Premise 3 is satisfied by construction.
    print("Premise 3 holds by construction of x1.")

    # 6. Evaluate the conclusion (Option D): Total Variation is 0.
    # TV(x0) = sum |x0_u - x0_v| over all edges
    print("\n--- Evaluating Conclusion D ---")
    tv = 0
    tv_equation_str = "TV(x0) = "
    for i, (u, v) in enumerate(edges):
        diff = abs(x0[u] - x0[v])
        tv += diff
        tv_equation_str += f"|{x0[u]} - {x0[v]}|"
        if i < len(edges) - 1:
            tv_equation_str += " + "
    
    tv_calc_str = "       = "
    for i, (u, v) in enumerate(edges):
        diff = abs(x0[u] - x0[v])
        tv_calc_str += str(diff)
        if i < len(edges) - 1:
            tv_calc_str += " + "
    
    print(tv_equation_str)
    print(tv_calc_str)
    print(f"       = {tv}")
    print("Conclusion: The total variation of the signal x0 on graph G is 0. This matches option D.")


if __name__ == '__main__':
    demonstrate_inference()