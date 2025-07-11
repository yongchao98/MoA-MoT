def solve():
    """
    This function determines the complexity categories for Luby's MIS algorithm on three classes of graphs.
    """
    # Analysis for f_1(n) on a cycle of length n
    # The graph is 2-regular. The complexity is known to be O(log n) and Omega(log n).
    # So, f_1(n) is Theta(log n). This falls into category 9.
    d1 = 9

    # Analysis for f_2(n) on any tree on n vertices of degree at most 100
    # The maximum degree is bounded by a constant (100). The complexity is O(log n).
    # A path is a tree, and it provides an Omega(log n) lower bound.
    # So, f_2(n) is Theta(log n). This falls into category 9.
    d2 = 9

    # Analysis for f_3(n) on any graph on n vertices of degree at most 100
    # The maximum degree is bounded by a constant (100). The complexity is O(log n).
    # There exist bounded-degree graphs (e.g., expanders) that provide an Omega(log n) lower bound.
    # So, f_3(n) is Theta(log n). This falls into category 9.
    d3 = 9
    
    final_answer = f"{d1}{d2}{d3}"
    
    print(f"The complexity for a cycle is f_1(n) = Theta(log n), which is category {d1}.")
    print(f"The complexity for a tree with max degree 100 is f_2(n) = Theta(log n), which is category {d2}.")
    print(f"The complexity for a graph with max degree 100 is f_3(n) = Theta(log n), which is category {d3}.")
    print(f"The three digit code is: {final_answer}")

solve()
print("<<<999>>>")