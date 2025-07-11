import math

def solve_labeling_number():
    """
    Determines the global labeling number of the graph K_{1,100}.
    """
    # The graph K_{1,100} is a star graph with one central vertex 'c' and 100 leaf vertices l_i.
    # The edges are e_i = (c, l_i), and we assign a positive integer label w_i to each edge.

    # The condition for a global labeling is that for any orientation of the edges and any
    # adjacent pair of vertices, the sums of their incoming edge labels must be different.
    # For K_{1,100}, the only adjacent pairs are (c, l_i).

    # By considering an orientation where the edge e_i points from c to l_i, the incoming
    # sum at l_i is w_i. The incoming sum at c is the sum of labels w_j for all edges e_j
    # that are oriented towards c. Since we can choose the orientation of all other edges
    # (j != i) arbitrarily, the sum at c can be the sum of any subset of {w_j | j != i}.

    # This leads to the following condition on the set of labels W = {w_1, ..., w_100}:
    # For any label w_i in W, w_i cannot be equal to the sum of labels of any subset of W \ {w_i}.

    # To find the minimum possible maximum label (the global labeling number k), we construct
    # the set of labels greedily, ensuring w_1 < w_2 < ... < w_100.
    # This simplifies the condition: w_i cannot be a sum of a subset of {w_1, ..., w_{i-1}}.

    # The greedy construction that minimizes the labels is as follows:
    # w_1 = 1
    # w_2 = 2 (smallest integer > 1 that is not a sum of a subset of {1})
    # w_3 = 4 (smallest integer > 2 that is not a sum of a subset of {1, 2}, i.e., not 1, 2, or 3)
    # This reveals the pattern w_i = 2^(i-1).

    # For K_{1,100}, we have n = 100 labels.
    n = 100

    # The set of labels is {2^0, 2^1, 2^2, ..., 2^(n-1)}.
    # The global labeling number, k, is the largest label in this set, which is w_n.
    
    # The final equation is k = 2^(n-1).
    base = 2
    exponent = n - 1
    
    # Calculate the result
    result = base ** exponent
    
    print("The problem is to find the global labeling number of the complete bipartite graph K_1,100.")
    print("The analysis shows that the optimal set of 100 labels follows the pattern w_i = 2^(i-1).")
    print(f"The largest label, k, determines the global labeling number.")
    print("\n--- Final Equation ---")
    print(f"k = {base}^({n} - 1)")
    
    # As requested, outputting each number in the final equation.
    print(f"Base: {base}")
    print(f"Exponent: {exponent}")
    
    print("\n--- Result ---")
    print(f"The global labeling number of K_1,100 is: {result}")

solve_labeling_number()