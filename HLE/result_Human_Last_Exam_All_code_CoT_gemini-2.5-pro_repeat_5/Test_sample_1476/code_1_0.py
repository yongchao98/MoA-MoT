import numpy as np

def solve():
    """
    This function demonstrates the logic for solving the problem.
    It shows that for the given conditions to hold, the edge signal x1 must be 0,
    which in turn means the total variation is 0.
    """
    print("Step 1: Define the graph structure and signals.")
    # We use a square graph with 4 nodes and 4 edges, which contains a cycle.
    # Nodes: 0, 1, 2, 3
    # Edges: (0,1), (1,2), (2,3), (3,0)
    # We need to define the vertex-edge incidence matrix B1.
    # We'll use the orientation 0->1, 1->2, 2->3, 3->0.
    # B1 has shape (num_vertices, num_edges)
    # B1[v, e] = -1 if v is the tail of e, +1 if v is the head, 0 otherwise.
    B1 = np.array([
        [-1,  0,  0,  1],  # Node 0 is tail of e0, head of e3
        [ 1, -1,  0,  0],  # Node 1 is head of e0, tail of e1
        [ 0,  1, -1,  0],  # Node 2 is head of e1, tail of e2
        [ 0,  0,  1, -1]   # Node 3 is head of e2, tail of e3
    ])

    print("Step 2: Start with a non-constant vertex signal x0.")
    x0 = np.array([10, 20, 10, 20])
    print(f"Initial x0 = {x0}")

    print("\nStep 3: Define x1 based on the problem statement.")
    # x1_e = |x0_u - x0_v|
    # Edge 0 (0,1): |10-20|=10
    # Edge 1 (1,2): |20-10|=10
    # Edge 2 (2,3): |10-20|=10
    # Edge 3 (3,0): |20-10|=10
    x1 = np.array([10, 10, 10, 10])
    print(f"This gives x1 = {x1}")
    print("All elements of x1 are non-negative, as required.")

    print("\nStep 4: Check the two main conditions from the problem.")
    # Condition A: "No cycles with non-zero sum"
    # The sum of non-negative x1 values on the cycle (0,1,2,3) is 10+10+10+10 = 40.
    # This is non-zero, so this x1 violates Condition A.
    print("Condition A Check: 'No cycles with non-zero sum'")
    cycle_sum = np.sum(x1)
    print(f"Sum of x1 values on the cycle is {cycle_sum}.")
    if cycle_sum == 0:
        print(" -> Condition A is satisfied.")
    else:
        print(" -> Condition A is NOT satisfied.")

    # Condition B: B1 * x1 = 0
    divergence = B1 @ x1
    print("\nCondition B Check: B1 * x1 = 0")
    print(f"B1.dot(x1) = {divergence}")
    if np.all(divergence == 0):
        print(" -> Condition B is satisfied.")
    else:
        print(" -> Condition B is NOT satisfied.")

    print("\nStep 5: Find a signal that satisfies all conditions.")
    print("To satisfy Condition A, all x1 values on the cycle must be 0.")
    print("This means x1 must be the zero vector.")
    x1_final = np.array([0, 0, 0, 0])
    print(f"Let's set x1 = {x1_final}")
    
    print("\nRe-checking conditions with the new x1:")
    # Check Condition A again
    cycle_sum_final = np.sum(x1_final)
    print(f"Condition A sum is now {cycle_sum_final}, which is 0. Condition A is satisfied.")
    
    # Check Condition B again
    divergence_final = B1 @ x1_final
    print(f"Condition B check: B1.dot(x1) = {divergence_final}. Condition B is satisfied.")
    
    print("\nStep 6: Calculate the final inferred quantity.")
    # The only way to satisfy all conditions is if x1 is the zero vector.
    # This happens if x0 is constant, e.g., x0 = [10, 10, 10, 10].
    # Total Variation is the sum of all elements in x1.
    total_variation = np.sum(x1_final)
    
    print("The total variation is the sum of the absolute differences on the edges.")
    print("This is the sum of the elements in x1.")
    final_eq_str = " + ".join(map(str, x1_final))
    print(f"Total Variation = {final_eq_str} = {total_variation}")
    print("\nConclusion: The premises together imply that the total variation must be 0.")

solve()