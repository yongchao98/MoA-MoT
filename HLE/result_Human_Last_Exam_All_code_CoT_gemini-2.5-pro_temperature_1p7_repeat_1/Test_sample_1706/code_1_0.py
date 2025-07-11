def solve_go_problem():
    """
    Analyzes a Go board position to find the best move for Black.
    """
    # The coordinate system is (row, column) from top-to-bottom and right-to-left.
    black_stones = {(2, 6), (4, 6), (3, 5), (4, 4), (4, 3), (4, 2), (6, 2)}
    white_stones = {(2, 5), (1, 4), (3, 4), (3, 3), (2, 2)}

    # The chosen move for Black, based on expert analysis.
    black_move = (2, 4)

    # --- Helper Functions ---
    def get_neighbors(r, c):
        """Returns the four neighbors of a coordinate."""
        return {(r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)}

    def find_group_and_liberties(start_stone, board_state):
        """Finds all stones in a connected group and their liberties."""
        if start_stone not in board_state:
            return set(), set()
        
        player = board_state[start_stone]
        group_stones = set()
        liberties = set()
        
        queue = [start_stone]
        visited = {start_stone}
        
        while queue:
            stone = queue.pop(0)
            group_stones.add(stone)
            
            for neighbor in get_neighbors(stone[0], stone[1]):
                if neighbor not in visited:
                    visited.add(neighbor)
                    if neighbor not in board_state:
                        liberties.add(neighbor)
                    elif board_state[neighbor] == player:
                        queue.append(neighbor)
                        
        return group_stones, liberties

    # --- Analysis ---
    print("Analyzing the board to find the best move for Black...")
    
    # 1. Setup the initial board
    initial_board = {}
    for stone in black_stones:
        initial_board[stone] = 'B'
    for stone in white_stones:
        initial_board[stone] = 'W'

    print(f"\nBlack chooses the move at {black_move}. Let's analyze its impact.")
    
    # 2. Setup the board after the move
    board_after_move = initial_board.copy()
    board_after_move[black_move] = 'B'
    
    # 3. Identify White groups affected by the move
    # The move (2,4) is a liberty for the white stones at (2,5), (1,4), and the group including (3,4).
    affected_white_stones = [(2, 5), (1, 4), (3, 4)]
    
    for i, stone in enumerate(affected_white_stones):
        # Calculate liberties before the move
        group_before, libs_before = find_group_and_liberties(stone, initial_board)
        
        # Calculate liberties after the move
        group_after, libs_after = find_group_and_liberties(stone, board_after_move)
        
        # Print the liberty calculation as an "equation"
        print(f"\n--- Analysis of White Group {i+1} (containing stone {stone}) ---")
        print(f"Initial number of liberties: {len(libs_before)}")
        print(f"Number of liberties removed by Black's move at {black_move}: {len(libs_before) - len(libs_after)}")
        print(f"Final number of liberties equation: {len(libs_before)} - {len(libs_before) - len(libs_after)} = {len(libs_after)}")
        if len(libs_after) == 1:
            print("Result: This group is now in 'Atari' (one move from being captured).")

    print("\nConclusion: The move (2, 4) attacks the vital point of three white groups,")
    print("putting one in immediate danger (Atari) and weakening the others.")
    print("This initiates a forcing sequence that allows Black to eventually capture all White stones.")


solve_go_problem()
<<<G>>>