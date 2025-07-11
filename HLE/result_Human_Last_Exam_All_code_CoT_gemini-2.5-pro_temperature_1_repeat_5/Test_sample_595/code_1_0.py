import math

def solve():
    """
    Calculates the maximum number of coordinate grid squares a triangle's perimeter can pass through.
    """
    # A small irrational offset to ensure no lattice points are on the perimeter
    epsilon = 1e-9

    # For a leg of length 18, we have dx^2 + dy^2 = 18^2 = 324.
    # To maximize floor(dx) + floor(dy), we find the best integer parts.
    # By checking integers from 1 to 17 for floor(dx), we find the max
    # sum is 25, achieved with floor(dx)=14 and floor(dy)=11.
    # dx=14, dy^2 = 324 - 14^2 = 128. dy = sqrt(128).
    dx1 = 14.0
    dy1 = math.sqrt(128) # approx 11.31

    # Place the right-angle vertex V2 at (epsilon, epsilon)
    v2 = (epsilon, epsilon)
    # The other vertices are positioned based on the leg displacements.
    # V1 is such that V2-V1 = (dx1, dy1)
    v1 = (epsilon - dx1, epsilon - dy1)
    # V3 is such that V3-V2 is perpendicular to V2-V1 and has length 18.
    # So V3-V2 = (-dy1, dx1)
    v3 = (epsilon - dy1, epsilon + dx1)

    # Function to calculate squares crossed by a segment from p1 to p2
    def count_squares(p1, p2):
        x1, y1 = p1
        x2, y2 = p2
        x_min, x_max = min(x1, x2), max(x1, x2)
        y_min, y_max = min(y1, y2), max(y1, y2)
        
        # N = 1 + num_vertical_crossings + num_horizontal_crossings
        num_squares = 1 + (math.floor(x_max) - math.floor(x_min)) + \
                      (math.floor(y_max) - math.floor(y_min))
        return num_squares

    # Calculate squares for each side
    # Leg 1: V1 to V2
    n1 = count_squares(v1, v2)
    # Leg 2: V2 to V3
    n2 = count_squares(v2, v3)
    # Hypotenuse: V3 to V1
    n3 = count_squares(v3, v1)
    
    # Overlaps: The 3 squares containing the vertices are counted twice.
    # Since we chose epsilon to be a small irrational number, the vertices
    # are in the interior of squares, not on edges or corners.
    overlaps = 3
    
    # Total unique squares
    k = n1 + n2 + n3 - overlaps
    
    print("Calculation for the optimal configuration:")
    print(f"Squares crossed by leg 1 (V1V2): N1 = {n1}")
    print(f"Squares crossed by leg 2 (V2V3): N2 = {n2}")
    print(f"Squares crossed by hypotenuse (V3V1): N3 = {n3}")
    print(f"Overlapping squares at vertices: {overlaps}")
    print(f"Total unique squares k = N1 + N2 + N3 - Overlaps")
    print(f"k = {n1} + {n2} + {n3} - {overlaps} = {k}")

solve()
<<<80>>>