import math

def count_crossings(p1, p2):
    """
    Counts the number of vertical and horizontal grid lines crossed by a line segment.
    """
    x1, y1 = p1
    x2, y2 = p2
    
    # Integers strictly between the x-coordinates
    min_x, max_x = sorted([x1, x2])
    # The number of integers n with min_x < n < max_x
    nv = math.floor(max_x - 1e-9) - math.ceil(min_x + 1e-9) + 1
    
    # Integers strictly between the y-coordinates
    min_y, max_y = sorted([y1, y2])
    # The number of integers m with min_y < m < max_y
    nh = math.floor(max_y - 1e-9) - math.ceil(min_y + 1e-9) + 1
    
    return nv + nh

def solve_triangle_grid_problem():
    """
    Calculates the maximum number of grid squares the perimeter of the triangle can cross.
    """
    # Side lengths of the isosceles right triangle
    leg_length = 18.0
    
    # Optimal orientation found by maximizing the projection perimeter
    # tan(theta) = 1/2, so sin(theta) = 1/sqrt(5), cos(theta) = 2/sqrt(5)
    
    # Define vectors for the two legs from the right-angle vertex
    v_leg1 = (leg_length * 2 / math.sqrt(5), leg_length * 1 / math.sqrt(5))
    v_leg2 = (-leg_length * 1 / math.sqrt(5), leg_length * 2 / math.sqrt(5))

    # Place the vertices with a small offset to avoid lattice points
    offset = 0.1  # A small non-integer offset
    
    # Right-angle vertex A
    A = (offset, offset)
    
    # Other two vertices B and C
    B = (A[0] + v_leg1[0], A[1] + v_leg1[1])
    C = (A[0] + v_leg2[0], A[1] + v_leg2[1])

    # Calculate crossings for each side of the triangle
    crossings_ab = count_crossings(A, B)
    crossings_ac = count_crossings(A, C)
    crossings_bc = count_crossings(B, C)
    
    total_crossings = crossings_ab + crossings_ac + crossings_bc

    print(f"Optimal placement of the triangle with side lengths {leg_length}, {leg_length}, {leg_length}*sqrt(2):")
    print(f"Vertex A at ({A[0]:.3f}, {A[1]:.3f})")
    print(f"Vertex B at ({B[0]:.3f}, {B[1]:.3f})")
    print(f"Vertex C at ({C[0]:.3f}, {C[1]:.3f})")
    print("-" * 30)
    print(f"Crossings for leg AB = {crossings_ab}")
    print(f"Crossings for leg AC = {crossings_ac}")
    print(f"Crossings for hypotenuse BC = {crossings_bc}")
    print("-" * 30)
    print(f"Total number of squares crossed = {crossings_ab} + {crossings_ac} + {crossings_bc} = {total_crossings}")
    print("-" * 30)
    print(f"The largest number k is: {total_crossings}")

solve_triangle_grid_problem()
