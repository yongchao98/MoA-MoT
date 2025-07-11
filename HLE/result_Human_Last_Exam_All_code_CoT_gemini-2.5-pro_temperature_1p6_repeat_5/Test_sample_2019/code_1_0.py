import itertools

def count_non_attacking_queens():
    """
    Calculates the number of ways to place 4 non-attacking queens
    on a 5x5 toroidal chessboard.
    """
    board_size = 5
    num_queens = 4

    # Generate all possible squares on the board (0-indexed)
    squares = [(r, c) for r in range(board_size) for c in range(board_size)]

    # Generate all possible combinations of placing the queens
    possible_placements = itertools.combinations(squares, num_queens)

    valid_placement_count = 0

    # Iterate through each combination of queen placements
    for placement in possible_placements:
        is_valid = True
        # Check for attacks between any pair of queens in the current placement
        for queen1, queen2 in itertools.combinations(placement, 2):
            r1, c1 = queen1
            r2, c2 = queen2

            # Check for attacks based on toroidal geometry
            if (r1 == r2 or                                           # Same row
                c1 == c2 or                                           # Same column
                (r1 - c1) % board_size == (r2 - c2) % board_size or   # Same main diagonal
                (r1 + c1) % board_size == (r2 + c2) % board_size):   # Same anti-diagonal
                
                is_valid = False
                break  # This pair attacks, so the placement is invalid
        
        if is_valid:
            valid_placement_count += 1
    
    # The final equation includes the parameters of the problem and the result.
    print(f"Board size: {board_size}x{board_size}")
    print(f"Number of queens: {num_queens}")
    print(f"Number of ways to place the non-attacking queens is: {valid_placement_count}")

# Execute the function to find and print the solution
count_non_attacking_queens()
