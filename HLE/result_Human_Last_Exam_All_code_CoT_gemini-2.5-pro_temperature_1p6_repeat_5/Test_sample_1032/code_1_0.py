import math

def solve_horse_area():
    """
    Calculates the reachable area for the horse.
    """
    
    # 1. Initial parameters
    rope_length = 7.0 / 2.0
    
    # 2. Calculate the total unobstructed area
    # The area of the region |x| + |y| <= L is 2 * L^2
    total_area = 2 * (rope_length ** 2)
    
    print(f"Step 1: Calculate the maximum possible area without any obstacles.")
    print(f"The rope has a taxi-cab length L = 7/2 = {rope_length}.")
    print(f"The area of the diamond shape defined by |x| + |y| <= L is 2 * L^2.")
    print(f"Maximum Area = 2 * ({rope_length})^2 = {total_area}\n")

    # 3. Account for the house area
    # The house consists of three unit squares.
    house_area = 3.0
    
    print(f"Step 2: Subtract the area of the house.")
    print(f"The house is made of three unit squares, so its area is {house_area}.\n")

    # 4. Calculate the unreachable "shadow" area
    # A shadow is cast behind the inner corner (-1, -1). This shadow lies within
    # the 1x1 square with corners at (-2,-1), (-1,-1), (-1,-2), and (-2,-2).
    # A point P(x,y) in this square is unreachable if the shortest path to it is > 3.5.
    # The shortest path is min(path_via_F, path_via_D), where F=(-2,-1) and D=(-1,-2).
    # length(path_via_F) = 3 + |x+2| + |y+1|
    # length(path_via_D) = 3 + |x+1| + |y+2|
    # A point is unreachable if BOTH paths are > 3.5.
    # This leads to the condition: -0.5 < x-y < 0.5.
    # The reachable parts of this 1x1 square are two small triangles at the corners.
    
    shadow_square_side = 1.0
    shadow_square_area = shadow_square_side ** 2
    
    # The vertices of the reachable triangles within the shadow square are:
    # Triangle 1: (-2, -1), (-1.5, -1), (-2, -1.5)
    # Triangle 2: (-1, -2), (-1.5, -2), (-1, -1.5)
    # Both are right triangles with legs of length 0.5.
    triangle_leg = 0.5
    reachable_triangle_area = 0.5 * triangle_leg * triangle_leg
    
    shadow_area = shadow_square_area - 2 * reachable_triangle_area
    
    print(f"Step 3: Calculate the area of the unreachable shadow region.")
    print(f"A shadow is cast behind the house's inner corner.")
    print(f"This shadow is the area of a {shadow_square_side}x{shadow_square_side} square minus two small reachable triangles.")
    print(f"Each reachable triangle has an area of 0.5 * {triangle_leg} * {triangle_leg} = {reachable_triangle_area}.")
    print(f"Unreachable Shadow Area = {shadow_square_area} - 2 * {reachable_triangle_area} = {shadow_area}\n")

    # 5. Final calculation
    final_area = total_area - house_area - shadow_area
    
    print(f"Step 4: Calculate the final area.")
    print(f"The final area is the maximum area minus the house area and the shadow area.")
    print(f"Final Area = {total_area} - {house_area} - {shadow_area} = {final_area}")

solve_horse_area()