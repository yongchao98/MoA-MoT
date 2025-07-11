import numpy as np

def solve_alien_invasion():
    """
    Solves the alien invasion chess board problem.

    This function simulates the growth of an alien colony on an 8x8 board.
    It starts with a specific initial configuration of 8 captured squares
    and applies the growth rule iteratively until the colony can no longer expand.
    Finally, it prints the maximal size of the colony found.
    """
    # Initialize an 8x8 board with all squares vacant (0).
    board = np.zeros((8, 8), dtype=int)

    # Chess coordinates are converted to 0-indexed (row, col) tuples.
    # 'a1' is (7, 0), 'h8' is (0, 7).
    # Rows are 1-8 (mapped to 7-0), columns are a-h (mapped to 0-7).
    def chess_to_coords(s):
        col = ord(s[0]) - ord('a')
        row = 8 - int(s[1])
        return row, col

    # The 8 initial captured squares.
    # d5, e5 are fixed.
    # d2, e2, b4, g4, b7, g7 are chosen for optimal growth.
    initial_squares_chess = ['d5', 'e5', 'd2', 'e2', 'b4', 'g4', 'b7', 'g7']
    
    initial_coords = [chess_to_coords(s) for s in initial_squares_chess]

    for r, c in initial_coords:
        board[r, c] = 1

    # Simulation loop
    while True:
        newly_captured = []
        # Find all vacant squares that can be captured
        for r in range(8):
            for c in range(8):
                if board[r, c] == 0:  # If the square is vacant
                    neighbor_count = 0
                    # Check the 4 neighbors (up, down, left, right)
                    for dr, dc in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
                        nr, nc = r + dr, c + dc
                        if 0 <= nr < 8 and 0 <= nc < 8 and board[nr, nc] == 1:
                            neighbor_count += 1
                    
                    if neighbor_count >= 2:
                        newly_captured.append((r, c))

        if not newly_captured:
            # If no new squares were captured, the process stops
            break
        
        # Capture all the new squares
        for r, c in newly_captured:
            board[r, c] = 1

    # The maximal size K is the total number of captured squares
    K = np.sum(board)
    
    print(f"The maximal size of the aliens' colony (K) is the sum of all captured squares.")
    print(f"The final number of captured squares is: {K}")

solve_alien_invasion()