import itertools

def solve_4_queens_on_5x5_toroidal():
    """
    This script calculates the number of ways to place 4 non-attacking queens on a 5x5 toroidal chessboard.

    A placement is valid if no two queens share the same row, column, or any wrapped diagonal.
    """
    
    # Define the parameters of the problem
    board_size = 5
    num_queens = 4

    # Initialize a counter for valid placements
    valid_placements_count = 0

    # Create a list of all squares on the board, represented by numbers from 0 to 24.
    # This is efficient for the combination function.
    all_squares = range(board_size * board_size)

    # Generate all possible combinations of placing `num_queens` on the `all_squares`.
    # This gives us every possible set of 4 squares.
    for placement_indices in itertools.combinations(all_squares, num_queens):
        
        # For each combination, we check if it's a valid placement.
        # First, convert the square indices back to (row, column) coordinates.
        placement_coords = [(i // board_size, i % board_size) for i in placement_indices]

        # Extract the properties of each queen's position
        rows = [coord[0] for coord in placement_coords]
        cols = [coord[1] for coord in placement_coords]
        # For a toroidal board, diagonals are calculated with modulo arithmetic
        diag1 = [(coord[0] - coord[1]) % board_size for coord in placement_coords]
        diag2 = [(coord[0] + coord[1]) % board_size for coord in placement_coords]

        # A placement is valid if all queens are on distinct rows, columns, and diagonals.
        # We can check this by seeing if the number of unique values for each property
        # is equal to the number of queens.
        if (len(set(rows)) == num_queens and
            len(set(cols)) == num_queens and
            len(set(diag1)) == num_queens and
            len(set(diag2)) == num_queens):
            valid_placements_count += 1
            
    # The problem statement requires printing the numbers involved.
    print(f"Board size N = {board_size}")
    print(f"Number of queens K = {num_queens}")
    print(f"Result: There are {valid_placements_count} ways to place {num_queens} non-attacking queens on a {board_size}x{board_size} toroidal chessboard.")

# Execute the solver function
solve_4_queens_on_5x5_toroidal()