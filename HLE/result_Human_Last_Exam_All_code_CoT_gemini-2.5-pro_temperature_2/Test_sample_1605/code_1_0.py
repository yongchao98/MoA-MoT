def solve():
    """
    This function explains the reasoning and prints the final answer.
    
    Let D be the disconnection number. We are looking for the number of homeomorphism
    classes of compact connected metric spaces X for which D(X) = 4.

    The disconnection number D(X) is defined as m(X) + 1, where m(X) is the maximum
    size of a set of points S such that X \ S is connected.
    So, we need to find spaces where m(X) = 3.

    1. High-dimensional spaces (like spheres or disks) are too connected; their
       disconnection number is infinite. We should look at 1-dimensional spaces (graphs).

    2. Let G be a graph with k endpoints (vertices of degree 1). Removing these k points
       usually leaves the graph connected, so m(G) >= k, and D(G) >= k + 1.
       For D(G) = 4, we must have k <= 3.

    3. We examine graphs based on the number of endpoints, k:
       - k=0 or 1 (graphs with cycles): A simple circle has D=2. More complex graphs
         with cycles (like a figure-8 or a circle with a tail) have D > 4 because
         cycles provide redundant paths, making them harder to disconnect.
       - k=2 (a simple arc): D = 3.
       - k=3 (a simple triod, or Y-shape): Removing the 3 endpoints leaves the
         graph connected, so m(G) >= 3. Any 4 points chosen will disconnect the graph.
         Thus, m(G) = 3 and D(G) = 4. This is a valid class.
       - k>=4 (e.g., a 4-star or "X-shape"): m(G) >= k >= 4, so D(G) >= 5.

    4. The analysis suggests that only trees are likely candidates, and the only tree
       that satisfies the condition is the one with 3 endpoints, which is
       homeomorphic to the simple triod.

    Therefore, there is only one such homeomorphism class.
    """
    
    # The number of homeomorphism classes with disconnection number equal to four.
    number_of_classes = 1
    
    print(number_of_classes)

solve()