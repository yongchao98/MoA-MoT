import itertools

def solve_toroidal_queens():
    """
    Calculates the number of ways to place 4 non-attacking queens on a 5x5 toroidal chessboard.
    """
    board_dim = 5
    num_queens = 4
    board_size = board_dim * board_dim
    
    # Generate all possible squares on the board, represented by numbers 0-24.
    all_squares = range(board_size)
    
    # Generate all combinations of placing `num_queens` on the board.
    queen_placements_indices = itertools.combinations(all_squares, num_queens)
    
    valid_placements_count = 0
    
    for placement_indices in queen_placements_indices:
        # Convert indices (0-24) to (row, col) coordinates.
        queens = [(i // board_dim, i % board_dim) for i in placement_indices]
        
        is_valid = True
        # Check every pair of queens in the current placement for attacks.
        for q1, q2 in itertools.combinations(queens, 2):
            r1, c1 = q1
            r2, c2 = q2
            
            # Check for row, column, and the two types of toroidal diagonal attacks.
            if (r1 == r2 or
                c1 == c2 or
                (r1 - c1) % board_dim == (r2 - c2) % board_dim or
                (r1 + c1) % board_dim == (r2 + c2) % board_dim):
                
                is_valid = False
                break  # An attack was found, so this placement is invalid.
        
        if is_valid:
            valid_placements_count += 1
            
    print(f"Board size: {board_dim}x{board_dim}")
    print(f"Number of queens: {num_queens}")
    print(f"Total number of ways to place {num_queens} non-attacking queens on a {board_dim}x{board_dim} toroidal board is: {valid_placements_count}")

# Execute the function to find and print the solution.
solve_toroidal_queens()
<<<50>>>