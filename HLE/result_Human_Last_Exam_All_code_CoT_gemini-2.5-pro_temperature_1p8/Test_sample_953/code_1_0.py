def solve():
    """
    This function determines the complexity category for the given algorithm on three different classes of graphs.
    
    The analysis shows that for all three graph classes (cycles, bounded-degree trees, and bounded-degree graphs),
    the algorithm terminates in Theta(log n) steps with high probability.
    
    1. Cycle of length n: f_1(n) = Theta(log n) -> Category 9
    2. Tree with n vertices and degree at most 100: f_2(n) = Theta(log n) -> Category 9
    3. Graph with n vertices and degree at most 100: f_3(n) = Theta(log n) -> Category 9
    
    Category 9 is for f(n) = Omega(log n). Since we have an O(log n) upper bound as well,
    Theta(log n) is the tight bound, which is included in Omega(log n).
    """
    
    d1 = 9  # For cycles
    d2 = 9  # For trees with degree at most 100
    d3 = 9  # For graphs with degree at most 100
    
    # The question asks for three digits d1d2d3 that encode the range.
    # We will print the three digits together as a single number.
    print(f"{d1}{d2}{d3}")

solve()
<<<999>>>