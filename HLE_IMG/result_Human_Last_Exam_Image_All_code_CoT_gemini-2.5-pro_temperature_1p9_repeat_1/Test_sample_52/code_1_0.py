def solve_tiling_periods():
    """
    This function calculates and prints the tiling periods for the four cases
    based on geometric reasoning.
    """
    
    # The periods are determined by the symmetry of the shapes.
    # Period 1: Shape can tile by translation alone (parallelogram or centrally symmetric).
    # Period 2: Minimum period for an asymmetric shape.
    
    # Case 1: Triangle (13, 31, 23) -> Asymmetric -> Period = 2
    period1 = 2
    
    # Case 2: Quadrilateral (10, 4, 23, 31) -> Parallelogram -> Period = 1
    period2 = 1
    
    # Case 3: Irregular Polygon (5, 15, ..., 7) -> Asymmetric -> Period = 2
    period3 = 2
    
    # Case 4: 13-gon (4, 5, ..., 13) -> Asymmetric (odd vertices) -> Period = 2
    period4 = 2
    
    periods = [period1, period2, period3, period4]
    
    print(f"{periods[0]},{periods[1]},{periods[2]},{periods[3]}")

solve_tiling_periods()