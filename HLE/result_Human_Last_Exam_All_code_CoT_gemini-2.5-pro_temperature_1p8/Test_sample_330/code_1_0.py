import numpy as np

def solve_square_vertices():
    """
    Calculates the vertices of a square given one point on each side.
    The method uses complex numbers.
    """
    # Step 1: Represent the points as complex numbers.
    # The problem assumes the points are given in order around the square perimeter.
    p1 = complex(0.3511, 0.2027)
    p2 = complex(0.6753, 0.8303)
    p3 = complex(-0.2845, 0.9905)
    p4 = complex(-0.128, 0.2218)

    # Step 2: Calculate the center of one of the possible squares.
    # We choose the second formula for the center.
    m_numerator = p1 + 1j*p2 - p3 - 1j*p4
    m_denominator = 2 * (1 + 1j)
    m = m_numerator / m_denominator

    # Step 3: Calculate the average of the points.
    p_avg = (p1 + p2 + p3 + p4) / 4

    # The vector from the center to a reference point is z.
    z = p_avg - m

    # Step 4: Calculate the four vertices of the square.
    # The vertices are m+z, m-z, m+iz, m-iz.
    # One of these, m+z, is simply p_avg.
    # The other is 2m - p_avg.
    v1 = m + z  # This is equal to p_avg
    v2 = m + 1j*z
    v3 = m - z
    v4 = m - 1j*z
    
    vertices = [(v.real, v.imag) for v in [v1, v2, v3, v4]]
    
    # Step 5: Sort the vertices by their x-coordinate.
    vertices.sort(key=lambda coord: coord[0])
    
    # Step 6: Print the sorted vertices formatted to two decimal places.
    for x, y in vertices:
        print(f"({x:.2f},{y:.2f})")

solve_square_vertices()