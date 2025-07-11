def solve_puzzle():
    """
    Calculates the minimum number of moves to move the red piece to the top-left corner.
    """
    # Grid and piece properties
    grid_size = 4
    piece_rows = 2
    piece_cols = 2

    # Initial and target positions of the top-left corner of the piece
    initial_pos = (2, 2)
    target_pos = (0, 0)

    # Calculate vertical and horizontal displacement needed
    vertical_displacement = abs(target_pos[0] - initial_pos[0])
    horizontal_displacement = abs(target_pos[1] - initial_pos[1])

    # Calculate moves for vertical shift
    # We need to shift 'piece_cols' number of columns, each by 'vertical_displacement' steps.
    vertical_moves = piece_cols * vertical_displacement
    
    # Calculate moves for horizontal shift
    # We need to shift 'piece_rows' number of rows, each by 'horizontal_displacement' steps.
    horizontal_moves = piece_rows * horizontal_displacement

    # Total minimum moves is the sum of vertical and horizontal moves
    total_moves = vertical_moves + horizontal_moves
    
    # Print the explanation of the calculation
    print("This puzzle can be solved by shifting entire rows and columns.")
    print("The goal is to move the 2x2 red piece from the bottom-right to the top-left.")
    print("\nStep 1: Calculate the moves needed to shift the piece vertically.")
    print(f"The piece needs to move up by {vertical_displacement} units.")
    print(f"The piece occupies {piece_cols} columns, so we shift these {piece_cols} columns up {vertical_displacement} times each.")
    print(f"Vertical Moves = {piece_cols} columns * {vertical_displacement} shifts = {vertical_moves} moves.")

    print("\nStep 2: Calculate the moves needed to shift the piece horizontally.")
    print(f"The piece needs to move left by {horizontal_displacement} units.")
    print(f"The piece occupies {piece_rows} rows, so we shift these {piece_rows} rows left {horizontal_displacement} times each.")
    print(f"Horizontal Moves = {piece_rows} rows * {horizontal_displacement} shifts = {horizontal_moves} moves.")

    print("\nStep 3: Calculate the total minimum moves.")
    print(f"Total Moves = Vertical Moves + Horizontal Moves")
    print(f"Total Moves = {vertical_moves} + {horizontal_moves} = {total_moves}")

solve_puzzle()