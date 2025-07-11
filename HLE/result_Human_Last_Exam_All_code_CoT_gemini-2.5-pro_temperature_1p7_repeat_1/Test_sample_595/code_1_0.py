import math

def solve():
    """
    Calculates the maximum number of grid squares a triangle with sides 18, 18, 18*sqrt(2)
    can pass through without its perimeter containing any lattice points.
    """

    # A square is identified by its bottom-left integer coordinates (floor(x), floor(y))
    def get_crossed_squares(p1, p2):
        """
        Calculates the set of unique grid squares crossed by a line segment from p1 to p2.
        A square is represented by a tuple (ix, iy) of its bottom-left corner.
        """
        x1, y1 = p1
        x2, y2 = p2
        
        squares = set()
        
        dx = x2 - x1
        dy = y2 - y1

        # The number of steps for simulation should be large enough to hit every square.
        # This is proportional to the number of grid lines crossed.
        steps = int(abs(dx) + abs(dy)) * 4 + 2
        
        if steps == 2: # Handles zero-length dx, dy case
            squares.add((math.floor(x1), math.floor(y1)))
            return squares

        for i in range(steps + 1):
            t = i / steps
            x = x1 + t * dx
            y = y1 + t * dy
            squares.add((math.floor(x), math.floor(y)))
            
        return squares

    # A small non-integer offset to avoid vertices and edges on grid lines.
    # We choose it such that 2*EPS is also not an integer.
    EPS = 0.00012345

    a = 18

    # Case 1: Axis-aligned legs
    A1 = (EPS, EPS)
    B1 = (a + EPS, EPS)
    C1 = (EPS, a + EPS)

    s_ab1 = get_crossed_squares(A1, B1)
    s_ac1 = get_crossed_squares(A1, C1)
    s_bc1 = get_crossed_squares(B1, C1)
    total_squares1 = len(s_ab1.union(s_ac1, s_bc1))

    # Case 2: Legs at 45 degrees to axes
    d = a / math.sqrt(2)
    A2 = (d + EPS, EPS)
    B2 = (A2[0] + d, A2[1] + d)
    C2 = (A2[0] - d, A2[1] + d)

    s_ab2 = get_crossed_squares(A2, B2)
    s_ac2 = get_crossed_squares(A2, C2)
    s_bc2 = get_crossed_squares(B2, C2)
    total_squares2 = len(s_ab2.union(s_ac2, s_bc2))

    # Case 3: Optimal angle rotation (tan(theta) = 1/2)
    # Leg vector components v1=(u,v), v2=(-v,u)
    v = a / math.sqrt(5)
    u = 2 * a / math.sqrt(5)

    # Place vertex A to keep all coordinates positive
    A3 = (v + EPS, EPS)
    B3 = (A3[0] + u, A3[1] + v)
    C3 = (A3[0] - v, A3[1] + u)

    s_ab3 = get_crossed_squares(A3, B3)
    s_ac3 = get_crossed_squares(A3, C3)
    s_bc3 = get_crossed_squares(B3, C3)
    total_squares3 = len(s_ab3.union(s_ac3, s_bc3))
    
    k = max(total_squares1, total_squares2, total_squares3)
    print(f"The largest number of squares k is {k}.")

solve()