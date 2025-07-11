import math

def solve_emperor_tomb_problem():
    """
    Solves the optimization problem for Emperor Qin's artifacts.
    """
    # Step 1: Define all constants from the problem description.
    rect_w = 140  # cm
    rect_h = 110  # cm

    # A circle of 20cm radius is cut from a 40x40cm square.
    circle_footprint = 40
    square_footprint = 10

    # Characters per artifact based on problem interpretation.
    # Assumption: 1 character = 9 bagua symbols. A circle holds 999 symbols.
    chars_per_circle = 999 // 9
    chars_per_square = 4

    def calculate_best_packing(container_w, container_h):
        """
        Calculates the number of circles (M) and squares (N) that can be
        packed using a greedy "circles-first" strategy.
        """
        # First, pack as many circles (40x40 squares) as possible in a grid.
        m_cols = container_w // circle_footprint
        m_rows = container_h // circle_footprint
        num_circles = m_cols * m_rows

        occupied_w = m_cols * circle_footprint
        occupied_h = m_rows * circle_footprint

        # Then, pack squares (10x10) in the remaining L-shaped area.
        # This area is split into two rectangles for easier calculation.
        
        # Rectangle 1: Full width of the container, remaining height
        rem_h_area = container_h - occupied_h
        n1 = (container_w // square_footprint) * (rem_h_area // square_footprint)

        # Rectangle 2: Occupied height, remaining width of the container
        rem_w_area = container_w - occupied_w
        n2 = (rem_w_area // square_footprint) * (occupied_h // square_footprint)
        
        num_squares = n1 + n2
        
        return num_squares, num_circles

    # Step 2: Calculate for both material orientations.
    # Orientation 1: 140x110
    n1, m1 = calculate_best_packing(rect_w, rect_h)
    k1 = (n1 * chars_per_square) + (m1 * chars_per_circle)

    # Orientation 2: 110x140 (rotated)
    n2, m2 = calculate_best_packing(rect_h, rect_w)
    k2 = (n2 * chars_per_square) + (m2 * chars_per_circle)

    # Step 3: Determine the best outcome and print the results.
    if k1 >= k2:
        best_n, best_m, best_k = n1, m1, k1
    else:
        best_n, best_m, best_k = n2, m2, k2

    print("Maximizing the number of engraved Chinese characters.")
    print(f"The optimal number of squares (N) is: {best_n}")
    print(f"The optimal number of circles (M) is: {best_m}")
    print("\nThe final calculation for the total number of characters (K) is:")
    print(f"K = ({best_n} squares × {chars_per_square} chars/square) + ({best_m} circles × {chars_per_circle} chars/circle)")
    print(f"K = ({best_n * chars_per_square}) + ({best_m * chars_per_circle}) = {best_k}")
    
    print("\nThe final answer in the format N:M:K is:")
    print(f"{best_n}:{best_m}:{best_k}")

solve_emperor_tomb_problem()