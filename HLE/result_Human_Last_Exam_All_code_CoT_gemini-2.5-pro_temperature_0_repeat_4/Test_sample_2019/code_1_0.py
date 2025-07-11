import itertools

def solve_toroidal_queens():
    """
    Calculates the number of ways to place a given number of non-attacking queens
    on a toroidal (wrap-around) chessboard.
    """
    board_size = 5
    num_queens = 4
    valid_placements_count = 0

    # Create a list of all possible squares (coordinates) on the board
    all_squares = []
    for r in range(board_size):
        for c in range(board_size):
            all_squares.append((r, c))

    # Generate all combinations of placing 'num_queens' on the available squares
    possible_placements = itertools.combinations(all_squares, num_queens)

    # Iterate through each possible placement to check for validity
    for placement in possible_placements:
        is_valid = True
        # Check every pair of queens in the current placement
        for i in range(num_queens):
            for j in range(i + 1, num_queens):
                r1, c1 = placement[i]
                r2, c2 = placement[j]

                # Check for any type of attack
                # 1. Row or Column attack (should not happen with this method, but good practice)
                if r1 == r2 or c1 == c2:
                    is_valid = False
                    break
                # 2. Main diagonal attack (toroidal)
                if (r1 - c1) % board_size == (r2 - c2) % board_size:
                    is_valid = False
                    break
                # 3. Anti-diagonal attack (toroidal)
                if (r1 + c1) % board_size == (r2 + c2) % board_size:
                    is_valid = False
                    break
            
            if not is_valid:
                break
        
        if is_valid:
            valid_placements_count += 1

    # Print the final result in a descriptive sentence
    print(f"On a {board_size}x{board_size} toroidal board, there are {valid_placements_count} ways to place {num_queens} non-attacking queens.")

# Execute the function to find the solution
solve_toroidal_queens()