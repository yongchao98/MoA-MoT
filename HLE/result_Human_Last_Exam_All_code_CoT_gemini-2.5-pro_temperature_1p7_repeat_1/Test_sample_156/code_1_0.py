def solve():
    """
    Calculates the maximal size of the alien colony.

    The optimal strategy is to use the 8 initial squares to form a barrier
    that encloses the largest possible area on the board. The aliens will
    then capture all squares in this enclosed region.

    An L-shaped barrier is most efficient for enclosing a rectangular corner.
    A barrier of k+l squares can enclose an area of k*l. We have 8 squares,
    so we need to maximize k*l subject to k+l=8.

    The fixed squares d5=(3,3) and e5=(3,4) must be in the barrier. This
    constrains our choices for k and l. Analysis shows the best we can do is
    wall off a 3x5 rectangular area.

    The enclosed area is therefore 3 * 5 = 15 squares.
    The maximal colony size K is the sum of the initial 8 squares plus the
    15 squares from the enclosed area.
    """
    
    board_size = 8
    initial_squares_count = 8
    
    # After analyzing the constraints, the largest area that can be enclosed
    # by an 8-square barrier containing d5 and e5 is a 3x5 rectangle.
    enclosed_area_rows = 3
    enclosed_area_cols = 5
    
    enclosed_area = enclosed_area_rows * enclosed_area_cols
    
    # The maximal size K is the sum of initial squares and the enclosed area.
    K = initial_squares_count + enclosed_area

    print(f"The problem is to find the maximum size of the alien colony, K.")
    print(f"The optimal strategy is to use the 8 initial squares to form a barrier enclosing the largest possible area.")
    print(f"K = (number of initial squares) + (size of enclosed area)")
    print(f"The number of initial squares is {initial_squares_count}.")
    print(f"The largest possible enclosed area, given the constraints, is a rectangle of {enclosed_area_rows} by {enclosed_area_cols} squares.")
    print(f"The size of this enclosed area is {enclosed_area_rows} * {enclosed_area_cols} = {enclosed_area}.")
    print("\nThe final calculation for K is:")
    print(f"{initial_squares_count} + ({enclosed_area_rows} * {enclosed_area_cols}) = {K}")

solve()
