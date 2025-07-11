import math

def count_crossings(p1, p2):
    """Calculates the number of grid lines a line segment crosses."""
    x1, y1 = p1
    x2, y2 = p2
    
    # Using the formula: |floor(x2) - floor(x1)| + |floor(y2) - floor(y1)|
    u1, v1 = math.floor(x1), math.floor(y1)
    u2, v2 = math.floor(x2), math.floor(y2)
    
    crossings = abs(u2 - u1) + abs(v2 - v1)
    return crossings, u1, v1, u2, v2

def main():
    """
    Calculates the maximum number of squares the triangle's perimeter can pass through.
    """
    
    # --- Case 1: Legs parallel to coordinate axes ---
    # To avoid lattice points, we place the right-angle vertex B at non-integer coordinates.
    # Let B = (0.5, 0.5). For y = -x + c, we get c=1. This line has integer points.
    # To avoid this, choose coordinates whose sum is not an integer.
    # Let B = (0.5, 0.25). A = (0.5, 18.25), C = (18.5, 0.25).
    # Hypotenuse AC is on y - 0.25 = -1 * (x - 18.5) => y = -x + 18.75.
    # y + x = 18.75 (not an integer), so no lattice points on the hypotenuse.
    # The other sides are x=0.5 and y=0.25, which also don't have lattice points.
    
    B1 = (0.5, 0.25)
    A1 = (0.5, 18.25)
    C1 = (18.5, 0.25)

    k_A1B1, u_A1, v_A1, u_B1, v_B1 = count_crossings(A1, B1)[0:5]
    k_B1C1, u_B1, v_B1, u_C1, v_C1 = count_crossings(B1, C1)[0:5]
    k_C1A1, u_C1, v_C1, u_A1, v_A1 = count_crossings(C1, A1)[0:5]
    k1 = k_A1B1 + k_B1C1 + k_C1A1
    
    # --- Case 2: Legs at 45 degrees to axes (hypotenuse parallel to an axis) ---
    # The change in x and y for legs of length 18 is 18/sqrt(2) = 9*sqrt(2)
    s = 9 * math.sqrt(2)
    
    # We choose coordinates for B to maximize floor values and avoid lattice points.
    # Slopes are 1 and -1. Equations: y-y0 = x-x0 and y-y0 = -(x-x0).
    # This means y-x = y0-x0 and y+x = y0+x0.
    # We must choose (x0, y0) such that y0-x0 and y0+x0 are not integers.
    # Let's try to maximize the result. The total crossings k is given by
    # k = 2 * (floor(x0+s)-floor(x0-s) + floor(y0+s))
    # To maximize floor(y0+s), we want y0 to be close to 1.
    # To have floor(x0+s)-floor(x0-s) be ceil(2s), we choose x0=0.5.
    # Let's choose x0 = 0.5, y0 = 0.6.
    # y0-x0 = 0.1 (not integer), y0+x0 = 1.1 (not integer). Safe.
    
    x0, y0 = 0.5, 0.6
    
    B2 = (x0, y0)
    A2 = (x0 + s, y0 + s)
    C2 = (x0 - s, y0 + s)

    k_A2B2, u_A2, v_A2, u_B2, v_B2 = count_crossings(A2, B2)[0:5]
    k_B2C2, u_B2, v_B2, u_C2, v_C2 = count_crossings(B2, C2)[0:5]
    k_C2A2, u_C2, v_C2, u_A2, v_A2 = count_crossings(C2, A2)[0:5]
    k2 = k_A2B2 + k_B2C2 + k_C2A2

    print("Investigating two main orientations to find the maximum number of squares.\n")
    print(f"Case 1: Triangle legs are parallel to the axes.")
    print(f"Calculation gives k = {k1}.\n")
    
    print(f"Case 2: Triangle legs are at 45 degrees to the axes.")
    print(f"This orientation yields a higher value of k = {k2}.")
    print("\nHere is the detailed calculation for the optimal case (Case 2):")
    
    print(f"\nVertices (float): A={A2}, B={B2}, C={C2}")
    
    # Integer parts (floors) of coordinates
    print(f"Integer part of vertices: A=({u_A2}, {v_A2}), B=({u_B2}, {v_B2}), C=({u_C2}, {v_C2})")

    print("\nNumber of squares crossed by each side:")
    print(f"Side AB crosses: |{u_B2} - {u_A2}| + |{v_B2} - {v_A2}| = {abs(u_B2 - u_A2)} + {abs(v_B2 - v_A2)} = {k_A2B2}")
    print(f"Side BC crosses: |{u_C2} - {u_B2}| + |{v_C2} - {v_B2}| = {abs(u_C2 - u_B2)} + {abs(v_C2 - v_B2)} = {k_B2C2}")
    print(f"Side CA crosses: |{u_A2} - {u_C2}| + |{v_A2} - {v_C2}| = {abs(u_A2 - u_C2)} + {abs(v_A2 - v_C2)} = {k_C2A2}")
    
    print("\nTotal number of squares crossed is the sum of the crossings for each side.")
    print(f"Total k = {k_A2B2} + {k_B2C2} + {k_C2A2} = {k2}")

    print(f"\nThe largest number k is {max(k1, k2)}.")


if __name__ == '__main__':
    main()
