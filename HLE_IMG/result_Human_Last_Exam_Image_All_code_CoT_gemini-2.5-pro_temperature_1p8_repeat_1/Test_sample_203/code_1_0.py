import math

def solve_heesch_number():
    """
    Determines the Heesch number for the given polygons.
    
    1.  The Heesch number of a shape is the maximum number of times (k) it can
        be surrounded by layers of copies of itself.
    2.  If a shape can tile the plane, it can be surrounded by an infinite
        number of layers, so its Heesch number is defined as infinity.
    3.  The polygon in the figure is a hexagon with a pair of opposite sides
        (top and bottom) that are equal in length and parallel.
    4.  A mathematical theorem states that any hexagon with these properties
        can tile the plane. This is a classic result in tiling theory.
    5.  Therefore, the polygon shown can tile the plane. Its Heesch number
        is infinity.
    6.  Since all three polygons in the image are identical, they all share
        the same Heesch number.
    """
    
    # Using the infinity symbol for representation as requested.
    heesch_number = "∞"
    
    # There are three identical polygons.
    results = [heesch_number, heesch_number, heesch_number]
    
    # Print the answer in the required format.
    print(','.join(results))

solve_heesch_number()