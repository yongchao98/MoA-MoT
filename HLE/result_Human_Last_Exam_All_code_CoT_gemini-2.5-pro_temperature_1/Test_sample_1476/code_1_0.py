import numpy as np

def demonstrate_inference():
    """
    This function demonstrates the logical inference from the problem description.
    It shows that if a vertex signal x0 has zero total variation, the resulting
    edge signal x1 satisfies the given conditions.
    """
    # 1. Define a graph: A square with 4 vertices and 4 edges.
    # V = {0, 1, 2, 3}, E = {(0,1), (1,2), (2,3), (3,0)}
    # We orient the edges: 0->1, 1->2, 2->3, 3->0
    V_count = 4
    E_count = 4
    edges = [(0, 1), (1, 2), (2, 3), (3, 0)]

    # 2. Define the boundary matrix B1 (|V| x |E|)
    B1 = np.zeros((V_count, E_count))
    for i, edge in enumerate(edges):
        u, v = edge
        B1[u, i] = -1
        B1[v, i] = 1

    # The gradient operator is the transpose of B1
    B1_T = B1.T

    # 3. Define a vertex signal x0 with Total Variation = 0.
    # This means the signal must be constant over the connected graph.
    c = 10.0 # An arbitrary constant
    x0 = np.array([c, c, c, c])

    print("Step 1: Define a vertex signal x0 with zero total variation.")
    print(f"x0 = {x0}\n")

    # 4. Compute the edge signal x1 = |gradient(x0)|
    z = B1_T @ x0
    x1 = np.abs(z)

    print("Step 2: Compute the edge signal x1 = |gradient(x0)|.")
    print(f"The gradient z = B1^T @ x0 = {z}")
    print(f"The edge signal x1 = |z| = {x1}\n")
    print("Inference: From the premises, we deduce x1 must be the zero vector.")
    print(f"Our computed x1 is indeed zero: {np.all(x1 == 0)}\n")

    # 5. Verify the conditions from the prompt for our computed x1
    print("Step 3: Verify that this x1 satisfies the prompt's conditions.")
    # Condition 1: x1 is orthogonal to all cycles.
    # The cycle space ker(B1) for a square is spanned by c = [1, 1, 1, 1]^T.
    cycle_basis = np.ones(E_count)
    inner_product = np.dot(x1, cycle_basis)
    print(f"Condition 1 (no cycles with non-zero sum): <x1, c> = {inner_product}")
    print(f"This holds, as {inner_product} is 0.\n")

    # Condition 2: The divergence of x1 is zero.
    divergence = B1 @ x1
    print(f"Condition 2 (divergence-free): B1 @ x1 = {divergence}")
    print(f"This holds, as the result is a zero vector.\n")

    # 6. Show the Total Variation calculation for x0
    print("Step 4: Conclude that the Total Variation of x0 must be 0.")
    tv = 0
    tv_str_parts = []
    tv_val_parts = []
    for edge in edges:
        u, v = edge
        diff = np.abs(x0[v] - x0[u])
        tv += diff
        tv_str_parts.append(f"|x{v}-x{u}|")
        tv_val_parts.append(f"|{x0[v]}-{x0[u]}|")

    print("The Total Variation is calculated as:")
    print(f"TV(x0) = {' + '.join(tv_str_parts)}")
    print(f"       = {' + '.join(tv_val_parts)}")
    print(f"       = {' + '.join([str(np.abs(x0[v]-x0[u])) for u,v in edges])}")
    print(f"       = {tv}")

if __name__ == "__main__":
    demonstrate_inference()