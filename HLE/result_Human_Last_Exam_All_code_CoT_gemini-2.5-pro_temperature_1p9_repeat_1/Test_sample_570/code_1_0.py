import math

def solve():
    """
    Calculates the minimal area of a convex domain that intersects all lines px+qy=1,
    where p and q are coprime integers.
    
    The problem reduces to finding the area of the convex body K defined by |x| + |y| <= 1.
    This body is a square with vertices at (1, 0), (0, 1), (-1, 0), and (0, -1).
    We use the Shoelace formula to calculate the area of this polygon.
    """
    
    # Vertices of the polygon |x| + |y| <= 1
    vertices = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    
    # Apply the Shoelace formula
    n = len(vertices)
    area = 0.0
    
    print("The minimal area corresponds to the area of the polygon with vertices at (1,0), (0,1), (-1,0), and (0,-1).")
    print("We calculate this area using the Shoelace formula: Area = 0.5 * |sum(x[i]*y[i+1] - x[i+1]*y[i])|")
    
    sum_parts = []
    for i in range(n):
        j = (i + 1) % n
        x_i, y_i = vertices[i]
        x_j, y_j = vertices[j]
        
        term = x_i * y_j - x_j * y_i
        sum_parts.append(term)
        area += term
        
    final_area = 0.5 * abs(area)

    # Outputting each number in the final equation as requested.
    print(f"\nThe equation is: Area = 0.5 * |({sum_parts[0]}) + ({sum_parts[1]}) + ({sum_parts[2]}) + ({sum_parts[3]})|")
    print(f"Area = 0.5 * |{area}|")
    print(f"Final Area = {final_area}")

solve()