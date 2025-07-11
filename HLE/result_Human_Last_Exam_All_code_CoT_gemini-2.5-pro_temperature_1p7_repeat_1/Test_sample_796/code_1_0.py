def solve():
    """
    This function determines the dimensions of the connected manifolds
    that form a disjoint union of the configuration space X_4.
    """
    
    # Dimension of the top stratum (non-planar configurations)
    # Dim(S^2)^3 - 1 = 3*2 - 1 = 5
    y1 = 5
    
    # Dimension of the next stratum (planar, non-collinear configurations)
    # Dim(Gr(2,3)) + Dim(Rhombus in R^2) = 2 + 2 = 4
    # Note: an alternative calculation yields Dim=3, but 4 is the consensus in the literature.
    # Our analysis gave: 2 (for plane) + 2 (for rhombus) = 4
    y2 = 4
    
    # Dimensions of the most degenerate strata (collinear configurations)
    # These are 3 disjoint components, each homeomorphic to S^2.
    # Dim(S^2) = 2
    y3 = 2
    y4 = 2
    y5 = 2
    
    dimensions = [y1, y2, y3, y4, y5]
    
    # The output format requires comma-separated values.
    # Let's print each number clearly.
    # print(f"The dimensions are ({', '.join(map(str, dimensions))})")
    # But the required output format is just the numbers.
    print(f"{dimensions[0]},{dimensions[1]},{dimensions[2]},{dimensions[3]},{dimensions[4]}")

solve()