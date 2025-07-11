import itertools

def solve_toroidal_queens():
    """
    Calculates the number of ways to place 4 non-attacking queens on a 5x5 toroidal chessboard.
    """
    N = 5
    NUM_QUEENS = 4
    
    # Create a list of all squares on the board, represented as (row, col) tuples.
    board_squares = [(r, c) for r in range(N) for c in range(N)]
    
    solution_count = 0
    
    # Generate all possible combinations of placing NUM_QUEENS on the board.
    possible_placements = itertools.combinations(board_squares, NUM_QUEENS)
    
    # Iterate through each combination to check if it's a valid solution.
    for placement in possible_placements:
        # A placement is a tuple of queen positions, e.g., ((r1, c1), (r2, c2), ...)
        
        # Use sets to efficiently check for uniqueness of rows, columns, and diagonals.
        rows = set()
        cols = set()
        main_diagonals = set()  # Identified by (r - c) mod N
        anti_diagonals = set()  # Identified by (r + c) mod N
        
        is_valid_placement = True
        for r, c in placement:
            # Check for row conflict
            if r in rows:
                is_valid_placement = False
                break
            rows.add(r)
            
            # Check for column conflict
            if c in cols:
                is_valid_placement = False
                break
            cols.add(c)
            
            # Check for main diagonal conflict (toroidal)
            main_diag_index = (r - c) % N
            if main_diag_index in main_diagonals:
                is_valid_placement = False
                break
            main_diagonals.add(main_diag_index)
            
            # Check for anti-diagonal conflict (toroidal)
            anti_diag_index = (r + c) % N
            if anti_diag_index in anti_diagonals:
                is_valid_placement = False
                break
            anti_diagonals.add(anti_diag_index)
            
        if is_valid_placement:
            solution_count += 1
            
    # The final "equation" is the result of this computation.
    # We print the final count, which is the solution.
    # The numbers in this final output are the board size, number of queens, and the total ways.
    print(f"Number of ways to place {NUM_QUEENS} non-attacking queens on a {N}x{N} toroidal board is: {solution_count}")

# Execute the function
solve_toroidal_queens()
<<<250>>>