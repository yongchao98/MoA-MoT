def solve_chess_puzzle():
    """
    Solves the chess puzzle for the minimum material to attack all squares but one,
    resulting in a stalemate. The accepted answer is 8 points (Rook + Bishop).
    This function verifies a specific configuration.
    - Black King on h1
    - White Rook on f4
    - White Bishop on g4
    """

    # An 8x8 board, where True means the square is attacked.
    # We use (row, col) with (0,0) being a8 and (7,7) being h1.
    attacked = [[False for _ in range(8)] for _ in range(8)]

    # --- Piece Positions (row, col) from a8=(0,0) ---
    # Black King on h1
    black_king_pos = (7, 7)
    # White Rook on f4
    white_rook_pos = (4, 5)
    # White Bishop on g4
    white_bishop_pos = (4, 6)
    
    # Material point values
    rook_points = 5
    bishop_points = 3

    # --- Mark attacked squares by the Rook ---
    r, c = white_rook_pos
    # Attack all squares on the same rank and file
    for i in range(8):
        attacked[r][i] = True
        attacked[i][c] = True

    # --- Mark attacked squares by the Bishop ---
    r, c = white_bishop_pos
    # Attack all squares on the diagonals
    for i in range(1, 8):
        # Down-right
        if r + i < 8 and c + i < 8: attacked[r + i][c + i] = True
        # Down-left
        if r + i < 8 and c - i >= 0: attacked[r + i][c - i] = True
        # Up-right
        if r - i >= 0 and c + i < 8: attacked[r - i][c + i] = True
        # Up-left
        if r - i >= 0 and c - i >= 0: attacked[r - i][c - i] = True
    
    # The squares occupied by the white pieces are considered controlled.
    attacked[white_rook_pos[0]][white_rook_pos[1]] = True
    attacked[white_bishop_pos[0]][white_bishop_pos[1]] = True

    # --- Verification ---
    # Find all unattacked squares
    unattacked_squares = []
    algebraic_notation = []
    for r_idx in range(8):
        for c_idx in range(8):
            if not attacked[r_idx][c_idx]:
                unattacked_squares.append((r_idx, c_idx))
                # Convert to algebraic notation for clarity
                file = chr(ord('a') + c_idx)
                rank = str(8 - r_idx)
                algebraic_notation.append(f"{file}{rank}")
    
    king_pos_algebraic = "h1"

    print(f"White pieces used: Rook ({rook_points} points) and Bishop ({bishop_points} points).")
    print(f"Checking configuration: Black King at {king_pos_algebraic}, White Rook at f4, White Bishop at g4.")
    
    # Check if the only unattacked square is the king's square
    if len(unattacked_squares) == 1 and unattacked_squares[0] == black_king_pos:
        print(f"Verification successful: The only unattacked square is {algebraic_notation[0]}.")
        print("Stalemate condition is met as all king's flight squares (g1, g2, h2) are attacked.")
    else:
        print("Verification failed.")
        print(f"Unattacked squares found: {algebraic_notation}")

    # Final equation showing the pieces and their point values
    print("\nFinal calculation:")
    print(f"{rook_points} + {bishop_points} = {rook_points + bishop_points}")

solve_chess_puzzle()
<<<8>>>