import itertools

def solve_toroidal_queens():
    """
    Calculates the number of ways to place K non-attacking queens on a
    toroidal NxN chessboard.
    """
    N = 5  # The size of the board (NxN)
    K = 4  # The number of queens to place

    # A list of all possible squares on the board, e.g., [(0, 0), (0, 1), ...]
    all_squares = list(itertools.product(range(N), repeat=2))

    # This will store the final count of valid placements
    valid_placements_count = 0

    # Generate all possible ways to choose K squares out of the N*N available squares
    for queen_placement in itertools.combinations(all_squares, K):
        is_valid = True
        # For a given placement of K queens, we must check every pair of queens
        # to see if they attack each other.
        for q1, q2 in itertools.combinations(queen_placement, 2):
            r1, c1 = q1
            r2, c2 = q2

            # Check for attacks based on toroidal rules
            # 1. Same row or same column
            if r1 == r2 or c1 == c2:
                is_valid = False
                break
            # 2. Same "difference" diagonal (wrapped)
            if (r1 - c1) % N == (r2 - c2) % N:
                is_valid = False
                break
            # 3. Same "sum" diagonal (wrapped)
            if (r1 + c1) % N == (r2 + c2) % N:
                is_valid = False
                break
        
        # If the inner loop completed without finding any attacks,
        # then this placement is valid.
        if is_valid:
            valid_placements_count += 1
    
    print(f"The number of ways to place {K} non-attacking queens on a {N}x{N} toroidal chessboard is: {valid_placements_count}")

# Execute the function to find the answer
solve_toroidal_queens()