import math

def calculate_minimal_area():
    """
    Calculates the minimal area of a convex domain that intersects all lines
    px+qy=1 for coprime integers p and q.

    The minimal domain is a hexagon that can be decomposed into a square and two triangles.
    """

    # 1. Define the vertices of the component shapes
    square_v1 = (1, 0)
    square_v2 = (0, 1)

    triangle_v1 = (0, 1)
    triangle_v2 = (1, 1)
    triangle_v3 = (1, 0)

    # 2. Calculate the area of the central square
    # The side length of the square is the distance between (1,0) and (0,1)
    side_length_sq = (square_v1[0] - square_v2[0])**2 + (square_v1[1] - square_v2[1])**2
    area_square = side_length_sq # Area = side_length^2

    # 3. Calculate the area of one of the triangles
    # The triangle is a right triangle. The lengths of the two legs can be calculated.
    leg1_length = math.sqrt((triangle_v2[0] - triangle_v1[0])**2 + (triangle_v2[1] - triangle_v1[1])**2)
    leg2_length = math.sqrt((triangle_v3[0] - triangle_v2[0])**2 + (triangle_v3[1] - triangle_v2[1])**2)
    area_triangle = 0.5 * leg1_length * leg2_length

    # The hexagon consists of one square and two identical triangles.
    total_area = area_square + 2 * area_triangle
    
    # 4. Print the result in the desired equation format
    # The user wants to see each number in the final equation.
    # The decomposition is one square of area 2.0 and two triangles of area 0.5 each.
    area_part1 = area_square
    area_part2 = area_triangle
    area_part3 = area_triangle
    final_area = area_part1 + area_part2 + area_part3
    
    print(f"The minimal area is the sum of the areas of a central square and two triangles.")
    print(f"Final Equation: {area_part1} + {area_part2} + {area_part3} = {final_area}")
    print(f"The minimal area is: {final_area}")

calculate_minimal_area()