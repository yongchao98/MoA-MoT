import numpy as np

def demonstrate_inference():
    """
    Demonstrates the logic to solve the problem using a simple triangle graph.
    """
    # Let's model a simple graph G: a triangle with nodes 0, 1, 2.
    # The edges are e0={0,1}, e1={1,2}, e2={2,0}.
    # We orient them as 0->1, 1->2, 2->0.
    nodes = [0, 1, 2]
    edges = [(0, 1), (1, 2), (2, 0)]
    
    # B1 is the vertex-edge incidence matrix.
    # Rows correspond to nodes (0, 1, 2), columns to edges (e0, e1, e2).
    B1 = np.array([
        [-1,  0,  1],  # Node 0 is start of e0, end of e2
        [ 1, -1,  0],  # Node 1 is end of e0, start of e1
        [ 0,  1, -1]   # Node 2 is end of e1, start of e2
    ])
    
    print("Step-by-step reasoning:\n")
    print("1. Premise: 'No cycles with non-zero sum' implies x^1 is a conservative flow.")
    print("   Mathematically: x^1 is in the image of B1.T (the space of gradients).")
    
    print("\n2. Premise: 'B1 * x^1 * 1^T = 0' implies B1 * x^1 = 0.")
    print("   Mathematically: x^1 is a solenoidal flow (divergence-free).")
    print("   x^1 is in the kernel (null space) of B1.\n")
    
    print("3. Key Insight from Hodge Theory:")
    print("   The space of conservative flows (im(B1.T)) and the space of solenoidal flows (ker(B1)) are orthogonal.")
    print("   The only vector that is in both a subspace and its orthogonal complement is the zero vector.\n")
    
    print("4. Conclusion so far:")
    print("   From (1) and (2), we must conclude that x^1 = 0.\n")
    
    print("5. Incorporating the third premise: x^1_e = |x^0_u - x^0_v|")
    print("   If x^1 = 0, then for every edge e={u,v}, |x^0_u - x^0_v| = 0.")
    print("   This means x^0 must be constant across all connected nodes.\n")
    
    print("6. Final Inference:")
    print("   The Total Variation of the signal x^0 is TV(x^0) = sum(|x^0_u - x^0_v|).")
    print("   Since each term is 0, the total sum must be 0.\n")
    print("----------------------------------------------------\n")
    print("Demonstration with a signal x^0 that results in Total Variation = 0:\n")
    
    # Define a vertex signal x^0 that is constant.
    x0 = np.array([10.0, 10.0, 10.0])
    print(f"Let's choose a vertex signal x^0 = {x0}")
    
    # Calculate the corresponding edge signal x^1
    x1_components = []
    print("\nCalculating the edge signal x^1 using x^1_e = |x^0_u - x^0_v|:")
    for u, v in edges:
        val = abs(x0[u] - x0[v])
        x1_components.append(val)
        print(f"  x^1_{{{u},{v}}} = |{x0[u]} - {x0[v]}| = {val}")
    x1 = np.array(x1_components)
    print(f"Resulting edge signal x^1 = {x1}\n")

    # Check if this x1 satisfies the premises
    print("Verifying premises for this x^1:")
    # Check premise 2: B1 @ x1 == 0
    div_x1 = B1 @ x1
    print(f"  Is x^1 in ker(B1)? Check: B1 @ x^1 = {div_x1}. This is the zero vector, so Yes.")
    # Check premise 1: is x1 in im(B1.T)? Yes, the zero vector is in every subspace.
    print(f"  Is x^1 in im(B1.T)? Yes, the zero vector is always in the image space.\n")
    
    # Calculate and print the final equation for Total Variation
    print("Calculating the Total Variation for x^0:")
    tv_sum_str = " + ".join([f"|{x0[u]} - {x0[v]}|" for u,v in edges])
    tv_val_str = " + ".join([f"{abs(x0[u]-x0[v])}" for u,v in edges])
    total_variation = np.sum(x1)
    
    print(f"TV(x^0) = {tv_sum_str}")
    print(f"        = {tv_val_str}")
    print(f"        = {total_variation}")
    
    print("\nThe calculation shows the total variation is 0.")
    print("This confirms that if the premises hold, the signal's total variation must be 0.")

demonstrate_inference()