import math

def solve():
    """
    Calculates the maximum number of coordinate grid squares a triangle with sides 18, 18, 18*sqrt(2)
    can pass through without its perimeter containing any lattice points.
    """

    # Step 1: Define the vertices for optimal placement.
    # We place the right-angle vertex at A, with legs parallel to the axes.
    # A small offset (0.5) is used to ensure no vertex or perimeter segment lies on a lattice point or line.
    A = (0.5, 0.5)
    B = (18.5, 0.5)
    C = (0.5, 18.5)

    # Step 2: Define a function to count squares crossed by a line segment.
    # Formula for a segment from P1 to P2: |floor(x2)-floor(x1)| + |floor(y2)-floor(y1)| + 1
    def squares_crossed_by_segment(p1, p2):
        x1, y1 = p1
        x2, y2 = p2
        # Use math.floor to get the integer part of the coordinates
        dx = abs(math.floor(x2) - math.floor(x1))
        dy = abs(math.floor(y2) - math.floor(y1))
        return dx + dy + 1

    # Step 3: Calculate the number of squares crossed by each side of the triangle.
    # Side AB (leg)
    n_ab = squares_crossed_by_segment(A, B)
    # Side AC (leg)
    n_ac = squares_crossed_by_segment(C, A) # The order doesn't matter for this calculation
    # Side BC (hypotenuse)
    n_bc = squares_crossed_by_segment(B, C)

    print(f"Number of squares for side AB (length 18): {n_ab}")
    print(f"Number of squares for side AC (length 18): {n_ac}")
    print(f"Number of squares for side BC (length 18*sqrt(2)): {n_bc}")
    
    # Step 4: Calculate the total number of unique squares.
    # The squares containing the vertices (A, B, C) are each counted twice,
    # so we subtract 3 from the total sum.
    k = n_ab + n_ac + n_bc - 3
    
    print(f"\nThe total number of squares k is calculated by summing the squares for each side and subtracting the 3 corner squares that are double-counted.")
    print(f"k = {n_ab} + {n_ac} + {n_bc} - 3 = {k}")

solve()