def solve_engraving_problem():
    """
    Calculates the optimal number of squares (N) and circles (M)
    to maximize the total number of engraved characters (K).
    """
    # --- 1. Define constants ---
    container_w, container_h = 140, 110
    square_side = 10
    circle_radius = 20
    circle_bounding_box = circle_radius * 2

    chars_per_square = 4
    chars_per_circle = 999

    # --- 2. Maximize M (number of circles) ---
    # The contribution of circles (999) is much higher than squares (4),
    # so we prioritize fitting as many circles as possible.
    m_fit_w = container_w // circle_bounding_box
    m_fit_h = container_h // circle_bounding_box
    M = m_fit_w * m_fit_h

    # --- 3. Calculate remaining space for N (squares) ---
    # The grid of circles occupies a rectangular area.
    circles_used_w = m_fit_w * circle_bounding_box
    circles_used_h = m_fit_h * circle_bounding_box

    # The remaining L-shaped area is split into two rectangles.
    # Rectangle 1
    rem1_w = container_w - circles_used_w
    rem1_h = container_h
    n1 = (rem1_w // square_side) * (rem1_h // square_side)

    # Rectangle 2
    rem2_w = circles_used_w
    rem2_h = container_h - circles_used_h
    n2 = (rem2_w // square_side) * (rem2_h // square_side)

    N = n1 + n2

    # --- 4. Calculate K (total characters) ---
    K = (chars_per_square * N) + (chars_per_circle * M)

    # --- 5. Output the result in the specified format ---
    # The final equation is K = 4*N + 999*M.
    # The output shows the value of each number N, M, and K.
    print(f"{N}:{M}:{K}")

solve_engraving_problem()