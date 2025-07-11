# Define the size of the board and the number of queens.
BOARD_SIZE = 5
NUM_QUEENS = 4

# This variable will store the total count of valid solutions.
solution_count = 0

def is_safe(row, col, placements):
    """
    Checks if it's safe to place a queen at (row, col) given the existing
    placements on the toroidal board.

    Args:
        row: The row of the potential new queen.
        col: The column of the potential new queen.
        placements: A list of (r, c) tuples for already placed queens.

    Returns:
        True if the placement is safe, False otherwise.
    """
    for placed_row, placed_col in placements:
        # Check for same row or same column attacks
        if row == placed_row or col == placed_col:
            return False
            
        # Check for toroidal diagonal attacks using modular arithmetic
        if (row - col) % BOARD_SIZE == (placed_row - placed_col) % BOARD_SIZE:
            return False
        if (row + col) % BOARD_SIZE == (placed_row + placed_col) % BOARD_SIZE:
            return False
            
    return True

def solve_queens_placements(queens_to_place, start_index, placements):
    """
    A recursive function to find all valid placements for the queens.

    Args:
        queens_to_place: The number of queens still needing to be placed.
        start_index: The board square index (0 to 24) to start searching from.
                     This ensures each combination is found only once.
        placements: A list of (r, c) coordinates of queens placed so far.
    """
    global solution_count
    
    # Base case: If no more queens are left to place, we've found a valid solution.
    if queens_to_place == 0:
        solution_count += 1
        return

    # Iterate through all possible squares starting from start_index
    for i in range(start_index, BOARD_SIZE * BOARD_SIZE):
        row = i // BOARD_SIZE
        col = i % BOARD_SIZE
        
        # If placing a queen at (row, col) is safe
        if is_safe(row, col, placements):
            # Recurse to place the next queen.
            # We add the current placement and search from the next square (i + 1).
            new_placements = placements + [(row, col)]
            solve_queens_placements(queens_to_place - 1, i + 1, new_placements)

def main():
    """
    Main function to initiate the search and print the result.
    """
    # Start the backtracking search from the first square (index 0) with no queens placed.
    solve_queens_placements(NUM_QUEENS, 0, [])
    
    print(f"Number of ways to place {NUM_QUEENS} non-attacking queens on a {BOARD_SIZE}x{BOARD_SIZE} toroidal chessboard: {solution_count}")

if __name__ == "__main__":
    main()
