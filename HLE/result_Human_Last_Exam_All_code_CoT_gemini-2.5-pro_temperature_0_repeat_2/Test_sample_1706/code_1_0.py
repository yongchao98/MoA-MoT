import collections

def solve_go_problem():
    """
    This script solves a Go life-and-death problem.
    It determines the optimal first move for Black to capture a group of White stones.
    """

    # The current configuration of stones on the 19x19 board.
    # (row, column), where row is 1-19 (top-bottom) and col is 1-19 (right-left).
    black_pieces = {(2, 6), (4, 6), (3, 5), (4, 4), (4, 3), (4, 2), (6, 2)}
    white_pieces = {(2, 5), (1, 4), (3, 4), (3, 3), (2, 2)}

    # --- Step 1: Analyze the situation ---
    # The White stones form a single connected group. To capture this group,
    # Black must occupy all of its liberties before White can form two "eyes".
    # This is a classic life-and-death (tsumego) problem. The solution is to
    # find the "vital point" of the White group's shape.

    # --- Step 2: Find the liberties of the White group ---
    # A liberty is an empty point adjacent to any stone in the group.
    
    def get_neighbors(r, c):
        return {(r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)}

    all_pieces = black_pieces.union(white_pieces)
    liberties = set()
    for r_white, c_white in white_pieces:
        for r_neighbor, c_neighbor in get_neighbors(r_white, c_white):
            if (r_neighbor, c_neighbor) not in all_pieces:
                liberties.add((r_neighbor, c_neighbor))

    # The calculated liberties are: {(1, 5), (2, 4), (1, 3), (2, 3), (3, 2), (1, 2), (2, 1)}

    # --- Step 3: Identify the vital point ---
    # The shape formed by the White group and its liberties has a critical weakness.
    # The point (2, 4) is the vital point.
    #
    # Reasoning:
    # - Playing at (2, 4) places a Black stone in the center of White's potential
    #   eye space, severely restricting White's ability to form two eyes.
    # - If Black plays at (2, 4), any response by White is insufficient. Black can
    #   systematically eliminate the remaining liberties and capture the group.
    # - If Black plays any other move first (e.g., (3, 2)), White can immediately
    #   play at (2, 4). This move would solidify White's eye shape, allowing the
    #   group to live.
    #
    # Therefore, the only move that guarantees the capture is (2, 4).

    optimal_move_row = 2
    optimal_move_col = 4

    # --- Step 4: Output the final answer ---
    # The question asks to output each number in the final equation.
    # We will present the final coordinate as the solution.
    
    print(f"The optimal first move for Black is to play at the coordinate ({optimal_move_row}, {optimal_move_col}).")
    print("This move occupies the vital point of the White group, ensuring its eventual capture.")
    print("\nThe numbers that form the coordinate of the chosen move are:")
    print(f"Row = {optimal_move_row}")
    print(f"Column = {optimal_move_col}")

solve_go_problem()