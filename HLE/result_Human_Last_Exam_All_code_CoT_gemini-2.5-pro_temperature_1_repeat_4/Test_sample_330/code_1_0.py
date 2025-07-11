import cmath
import math

def solve_square_from_points():
    """
    Finds the vertices of a square given one point on each of its four sides.
    """
    # The four given points, assumed to be in cyclic order on the square's sides.
    # P1 on side AB, P2 on BC, P3 on CD, P4 on DA
    p1 = (0.3511, 0.2027)
    p2 = (0.6753, 0.8303)
    p3 = (-0.2845, 0.9905)
    p4 = (-0.128, 0.2218)
    
    points = [p1, p2, p3, p4]

    # Calculate differences needed for the quadratic equation coefficients
    # delta_x_13 = p1.x - p3.x
    delta_x_13 = p1[0] - p3[0]
    # delta_y_13 = p1.y - p3.y
    delta_y_13 = p1[1] - p3[1]
    # delta_x_24 = p2.x - p4.x
    delta_x_24 = p2[0] - p4[0]
    # delta_y_24 = p2.y - p4.y
    delta_y_24 = p2[1] - p4[1]

    # Coefficients of the quadratic equation Am^2 + Bm + C = 0 for the slope m
    A = delta_x_13**2 - delta_y_24**2
    B = -2 * (delta_x_13 * delta_y_13 + delta_x_24 * delta_y_24)
    C = delta_y_13**2 - delta_x_24**2

    print("The quadratic equation for the slope 'm' of the square's sides is:")
    print(f"{A:.4f} * m^2 + {B:.4f} * m + {C:.4f} = 0")
    print("")

    # Solve the quadratic equation for m
    discriminant = cmath.sqrt(B**2 - 4 * A * C)
    m1 = ((-B + discriminant) / (2 * A)).real
    m2 = ((-B - discriminant) / (2 * A)).real

    solutions = []
    
    for m in [m1, m2]:
        if abs(m) > 1e10: continue # Skip if slope is practically infinite
        
        # Calculate vertices for the current slope m
        # A = intersection of lines through P4 and P1
        # B = intersection of lines through P1 and P2
        # C = intersection of lines through P2 and P3
        # D = intersection of lines through P3 and P4
        
        # Denominator in vertex calculation formulas
        den = 1 + m**2
        
        # Vertex formulas derived from line intersections
        x_A = (m**2 * p1[0] - m * p1[1] + p4[0] + m * p4[1]) / den
        y_A = m * x_A - (m * p1[0] - p1[1])
        
        x_B = (m**2 * p1[0] - m * p1[1] + p2[0] + m * p2[1]) / den
        y_B = m * x_B - (m * p1[0] - p1[1])

        x_C = (m**2 * p3[0] - m * p3[1] + p2[0] + m * p2[1]) / den
        y_C = m * x_C - (m * p3[0] - p3[1])

        x_D = (m**2 * p3[0] - m * p3[1] + p4[0] + m * p4[1]) / den
        y_D = m * x_D - (m * p3[0] - p3[1])

        vertices = [(x_A, y_A), (x_B, y_B), (x_C, y_C), (x_D, y_D)]
        
        # Check if the given points lie on the segments of the square sides
        side_points = [p1, p2, p3, p4]
        side_vertices = [(vertices[0], vertices[1]), (vertices[1], vertices[2]), 
                         (vertices[2], vertices[3]), (vertices[3], vertices[0])]
        
        all_points_on_sides = True
        for i in range(4):
            p = side_points[i]
            v1 = side_vertices[i][0]
            v2 = side_vertices[i][1]
            
            # Check if point p's coordinates are between the vertex coordinates
            on_segment = (min(v1[0], v2[0]) - 1e-9 <= p[0] <= max(v1[0], v2[0]) + 1e-9 and
                          min(v1[1], v2[1]) - 1e-9 <= p[1] <= max(v1[1], v2[1]) + 1e-9)
            if not on_segment:
                all_points_on_sides = False
                break
        
        if all_points_on_sides:
            solutions.append(vertices)

    # There should be only one valid solution where points are on the segments
    if not solutions:
        print("No valid square found for the given point order.")
        return
        
    final_vertices = solutions[0]
    
    # Sort vertices by x-coordinate
    final_vertices.sort(key=lambda v: v[0])
    
    print("The coordinates of the vertices of the square are:")
    for v in final_vertices:
        print(f"({v[0]:.2f}, {v[1]:.2f})")
    
    # Format for the final answer block
    answer_string = ""
    for v in final_vertices:
        answer_string += f"({v[0]:.2f},{v[1]:.2f})"
    # The submission format seems to prefer just the values, let's provide just the text.
    answer_text_list = []
    for v in final_vertices:
        answer_text_list.append(f"({v[0]:.2f}, {v[1]:.2f})")
    
    final_answer = ", ".join(answer_text_list)
    print(f"\n<<<{final_answer}>>>")

solve_square_from_points()