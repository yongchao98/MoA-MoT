import math

def calculate_k():
    """
    Calculates the largest number k of coordinate grid squares the triangle's
    perimeter can pass through.
    """
    
    # Side lengths of the isosceles right triangle
    side = 18.0

    # Optimal orientation angle theta where tan(theta) = 1/2
    theta = math.atan(0.5)
    cos_theta = math.cos(theta)
    sin_theta = math.sin(theta)

    # The vectors for the two shorter sides from the right-angle vertex A
    # vector_AB = (side * cos_theta, side * sin_theta)
    # vector_AC = (-side * sin_theta, side * cos_theta)
    # Let A be at the origin (for simplicity, as integer part of coords cancels out)
    # We choose fractional parts of A's coordinate to be almost 1 to maximize k.
    
    # Coordinates of vertices B and C relative to A=(x0,y0)
    # where x0, y0 have fractional parts approaching 1.
    x_B_rel = side * cos_theta  # approx 16.099
    y_B_rel = side * sin_theta  # approx 8.049
    x_C_rel = -side * sin_theta # approx -8.049
    y_C_rel = side * cos_theta  # approx 16.099

    # Floored coordinates of the vertices A, B, C.
    # We assume A is at (ix + (1-eps), iy + (1-eps)) for integer ix, iy.
    # Then floor(xA) = ix, floor(yA) = iy.
    # floor(xB) = floor(ix + 1-eps + x_B_rel) = ix + floor(1 + x_B_rel) = ix + ceil(x_B_rel)
    # floor(yB) = iy + ceil(y_B_rel)
    # floor(xC) = floor(ix + 1-eps + x_C_rel) = ix + floor(x_C_rel)
    # floor(yC) = iy + ceil(y_C_rel)

    # Let integer parts of A's coordinates be 0 for simplicity.
    floor_xA = 0
    floor_yA = 0

    floor_xB = math.ceil(x_B_rel)
    floor_yB = math.ceil(y_B_rel)

    floor_xC = math.floor(x_C_rel)
    floor_yC = math.ceil(y_C_rel)

    # Calculate the Manhattan (L1) distances between the floored vertices
    l1_AB = abs(floor_xB - floor_xA) + abs(floor_yB - floor_yA)
    l1_AC = abs(floor_xC - floor_xA) + abs(floor_yC - floor_yA)
    l1_BC = abs(floor_xC - floor_xB) + abs(floor_yC - floor_yB)
    
    # The total number of squares k is the sum of these L1 distances
    k = l1_AB + l1_AC + l1_BC

    print(f"Calculation for the optimal orientation (tan(theta)=1/2):")
    print(f"Floored vertices relative to A's floored position:")
    print(f"A_floor = ({floor_xA}, {floor_yA})")
    print(f"B_floor = ({floor_xB}, {floor_yB})")
    print(f"C_floor = ({floor_xC}, {floor_yC})")
    print("-" * 20)
    print(f"Number of squares for side AB = |{floor_xB} - {floor_xA}| + |{floor_yB} - {floor_yA}| = {l1_AB}")
    print(f"Number of squares for side AC = |{floor_xC} - {floor_xA}| + |{floor_yC} - {floor_yA}| = {l1_AC}")
    print(f"Number of squares for side BC = |{floor_xC} - {floor_xB}| + |{floor_yC} - {floor_yB}| = {l1_BC}")
    print("-" * 20)
    print(f"Total number of squares k = {l1_AB} + {l1_AC} + {l1_BC} = {k}")

calculate_k()