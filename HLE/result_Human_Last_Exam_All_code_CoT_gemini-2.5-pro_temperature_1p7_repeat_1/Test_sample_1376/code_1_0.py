def solve_engraving_problem():
    """
    Calculates the optimal number of squares (N) and circles (M)
    to maximize the total number of engraved characters (K).
    """

    # 1. Define constants from the problem description
    material_w = 140  # cm
    material_h = 110  # cm
    square_side = 10  # cm
    circle_diameter = 40  # cm (20cm radius)
    chars_per_square = 4
    chars_per_circle = 999

    print("Objective: Maximize K = (4 * N) + (999 * M)")
    print("Strategy: Prioritize circles (M) as they provide more characters per area.\n")

    # 2. Calculate the maximum number of circles (M)
    # We find how many 40x40cm blocks fit along the material's width and height.
    m_along_w = material_w // circle_diameter
    m_along_h = material_h // circle_diameter
    M = m_along_w * m_along_h
    print(f"Calculating the maximum number of circles (M):")
    print(f"Circles along 140cm width = floor(140 / 40) = {m_along_w}")
    print(f"Circles along 110cm height = floor(110 / 40) = {m_along_h}")
    print(f"Total circles M = {m_along_w} * {m_along_h} = {M}\n")

    # 3. Calculate the number of squares (N) that fit in the remaining area
    # The M circles form a rectangular block.
    circles_block_w = m_along_w * circle_diameter
    circles_block_h = m_along_h * circle_diameter

    # The remaining area is an L-shape. We can divide it into two rectangles to calculate N.
    # Rectangle 1: The vertical part of the 'L'
    rem_rect1_w = material_w - circles_block_w
    rem_rect1_h = material_h
    n1 = (rem_rect1_w // square_side) * (rem_rect1_h // square_side)

    # Rectangle 2: The horizontal part of the 'L', not overlapping with Rect 1
    rem_rect2_w = circles_block_w
    rem_rect2_h = material_h - circles_block_h
    n2 = (rem_rect2_w // square_side) * (rem_rect2_h // square_side)

    N = n1 + n2
    print(f"Calculating the maximum number of squares (N) in the remaining space:")
    print(f"The {M} circles occupy a {circles_block_w}cm x {circles_block_h}cm area.")
    print(f"The remaining L-shape can be split into a {rem_rect1_w}x{rem_rect1_h}cm rectangle and a {rem_rect2_w}x{rem_rect2_h}cm rectangle.")
    print(f"Squares in first rectangle = ({rem_rect1_w}/{square_side}) * ({rem_rect1_h}/{square_side}) = {n1}")
    print(f"Squares in second rectangle = ({rem_rect2_w}/{square_side}) * ({rem_rect2_h}/{square_side}) = {n2}")
    print(f"Total squares N = {n1} + {n2} = {N}\n")

    # 4. Calculate the total characters (K)
    k_from_squares = N * chars_per_square
    k_from_circles = M * chars_per_circle
    K = k_from_squares + k_from_circles
    
    print("Calculating the total number of characters (K):")
    # Outputting each number in the final equation as requested
    print(f"Final Equation: {chars_per_square} * {N} + {chars_per_circle} * {M} = {k_from_squares} + {k_from_circles} = {K}\n")
    
    # 5. Print the final result in the specified N:M:K format
    print("Final Answer:")
    print(f"{N}:{M}:{K}")

solve_engraving_problem()