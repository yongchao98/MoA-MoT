import itertools

def solve_toroidal_n_queens():
    """
    Calculates the number of ways to place K non-attacking queens on an
    N x N toroidal chessboard by checking all possible placements.
    """
    N = 5  # The size of the board (NxN)
    K = 4  # The number of queens to place

    # Generate a list of all squares on the board, represented by indices 0 to (N*N - 1)
    all_squares_indices = range(N * N)

    # Function to convert a 1D index to 2D (row, col) coordinates
    def to_coords(index, board_size):
        row = index // board_size
        col = index % board_size
        return (row, col)

    # Function to check if a specific placement of queens is valid
    def is_valid_placement(placement, board_size):
        # Convert the chosen indices to (row, col) coordinates
        queen_coords = [to_coords(i, board_size) for i in placement]
        
        # Check all unique pairs of queens for attacks
        for q1, q2 in itertools.combinations(queen_coords, 2):
            r1, c1 = q1
            r2, c2 = q2
            
            # 1. Check for row attack
            if r1 == r2:
                return False
                
            # 2. Check for column attack
            if c1 == c2:
                return False
                
            # 3. Check for main diagonal attack (toroidal)
            if (r1 - c1) % board_size == (r2 - c2) % board_size:
                return False
            
            # 4. Check for anti-diagonal attack (toroidal)
            if (r1 + c1) % board_size == (r2 + c2) % board_size:
                return False
                
        # If no attacks are found, the placement is valid
        return True

    # Counter for valid placements
    valid_count = 0
    
    # Iterate through all combinations of K squares on the board
    for placement_indices in itertools.combinations(all_squares_indices, K):
        if is_valid_placement(placement_indices, N):
            valid_count += 1
            
    # As requested, output the numbers from the final result.
    # "The equation is: placing K queens on an NxN board results in X ways"
    print(f"Board Dimension (N) = {N}")
    print(f"Number of Queens (K) = {K}")
    print(f"Total Valid Placements = {valid_count}")

solve_toroidal_n_queens()