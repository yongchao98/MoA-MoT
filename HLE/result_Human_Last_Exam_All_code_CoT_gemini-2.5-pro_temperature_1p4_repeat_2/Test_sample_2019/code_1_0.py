import itertools

def solve_toroidal_queens():
    """
    Calculates the number of ways to place 4 non-attacking queens on a 5x5
    toroidal chessboard.
    """
    N = 5
    NUM_QUEENS = 4
    board_size = N * N
    
    # Generate all possible squares on the board, represented by indices 0-24.
    all_squares = range(board_size)
    
    # Generate all combinations of placing NUM_QUEENS queens on the board.
    possible_placements = itertools.combinations(all_squares, NUM_QUEENS)
    
    valid_solution_count = 0
    
    # Iterate through each possible placement of 4 queens.
    for placement in possible_placements:
        # Convert the 1D square indices to 2D (row, col) coordinates.
        coords = []
        for square_index in placement:
            row = square_index // N
            col = square_index % N
            coords.append((row, col))
            
        # Assume the current placement is valid until an attack is found.
        is_valid = True
        
        # Check all pairs of queens for attacks.
        for queen1, queen2 in itertools.combinations(coords, 2):
            r1, c1 = queen1
            r2, c2 = queen2
            
            # Check for attacks:
            # 1. Same row
            if r1 == r2:
                is_valid = False
                break
            # 2. Same column
            if c1 == c2:
                is_valid = False
                break
            # 3. Same main diagonal (toroidal)
            if (r1 - c1) % N == (r2 - c2) % N:
                is_valid = False
                break
            # 4. Same anti-diagonal (toroidal)
            if (r1 + c1) % N == (r2 + c2) % N:
                is_valid = False
                break
        
        # If the inner loop completed without finding an attack, it's a valid solution.
        if is_valid:
            valid_solution_count += 1
            
    # The problem asks to output the numbers in the final equation.
    # Here, we simply state the result of the calculation.
    print(f"The number of ways to place {NUM_QUEENS} non-attacking queens on a {N}x{N} toroidal chessboard is {valid_solution_count}.")

solve_toroidal_queens()