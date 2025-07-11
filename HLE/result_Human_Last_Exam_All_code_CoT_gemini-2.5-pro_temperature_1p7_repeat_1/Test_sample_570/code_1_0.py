import math

def solve():
    """
    Calculates the minimal area of a convex domain that intersects all lines
    px+qy=1 for coprime integers p, q.

    The reasoning shows that this minimal area is achieved by the convex set
    defined by |x| + |y| <= 1. This is a square rotated by 45 degrees,
    with vertices at (1,0), (0,1), (-1,0), and (0,-1).

    This function calculates the area of this polygon using the Shoelace formula.
    """
    # Vertices of the polygon
    vertices = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    
    # Extract x and y coordinates
    x = [v[0] for v in vertices]
    y = [v[1] for v in vertices]
    n = len(vertices)

    # Apply the Shoelace formula: Area = 0.5 * |(x1y2 + x2y3 + ... + xny1) - (y1x2 + y2x3 + ... + ynx1)|
    sum1 = 0
    sum2 = 0
    
    print("Applying the Shoelace formula for the polygon with vertices (1,0), (0,1), (-1,0), (0,-1):")
    print("Area = 0.5 * | (x1*y2 + x2*y3 + x3*y4 + x4*y1) - (y1*x2 + y2*x3 + y3*x4 + y4*x1) |")
    
    # Calculate sum1 and show its components
    print("\nCalculating the first part of the formula:")
    sum1_terms = []
    for i in range(n):
        j = (i + 1) % n
        term = x[i] * y[j]
        sum1 += term
        sum1_terms.append(f"{x[i]}*{y[j]}")
    print(f"  sum1 = {' + '.join(sum1_terms)}")
    print(f"  sum1 = {1*1} + {0*0} + {-1*(-1)} + {0*0} = {sum1}")
    
    # Calculate sum2 and show its components
    print("\nCalculating the second part of the formula:")
    sum2_terms = []
    for i in range(n):
        j = (i + 1) % n
        term = y[i] * x[j]
        sum2 += term
        sum2_terms.append(f"{y[i]}*{x[j]}")
    print(f"  sum2 = {' + '.join(sum2_terms)}")
    print(f"  sum2 = {0*0} + {1*(-1)} + {0*0} + {(-1)*1} = {sum2}")

    # Calculate the final area
    area = 0.5 * abs(sum1 - sum2)
    
    print("\nFinal Calculation:")
    print(f"Area = 0.5 * |{sum1} - ({sum2})|")
    print(f"Area = 0.5 * |{sum1 - sum2}|")
    print(f"Area = {area}")

solve()
<<<2>>>