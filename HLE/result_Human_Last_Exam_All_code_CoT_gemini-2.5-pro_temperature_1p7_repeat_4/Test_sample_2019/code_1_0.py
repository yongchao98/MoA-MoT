import itertools

def solve_toroidal_queens():
    """
    Calculates the number of ways to place a given number of non-attacking
    queens on a toroidal chessboard of a given size.
    """
    # --- Problem Parameters ---
    BOARD_SIZE = 5
    NUM_QUEENS = 4

    # --- Calculation ---

    # Create a list of all possible cell coordinates, e.g., (0, 0), (0, 1), ...
    cells = [(r, c) for r in range(BOARD_SIZE) for c in range(BOARD_SIZE)]

    # Generate all unique combinations of placing NUM_QUEENS on the board cells
    potential_placements = itertools.combinations(cells, NUM_QUEENS)

    solution_count = 0
    # Iterate through each combination of queen placements
    for placement in potential_placements:
        is_valid = True
        # Check all pairs of queens within the current placement for attacks
        for queen1, queen2 in itertools.combinations(placement, 2):
            r1, c1 = queen1
            r2, c2 = queen2

            # Check for attacks on the toroidal board
            # 1. Row attack
            if r1 == r2:
                is_valid = False
                break
            # 2. Column attack
            if c1 == c2:
                is_valid = False
                break
            # 3. Main diagonal attack (toroidal)
            if (r1 - c1) % BOARD_SIZE == (r2 - c2) % BOARD_SIZE:
                is_valid = False
                break
            # 4. Anti-diagonal attack (toroidal)
            if (r1 + c1) % BOARD_SIZE == (r2 + c2) % BOARD_SIZE:
                is_valid = False
                break
        
        # If no attacks were found after checking all pairs, it's a valid solution
        if is_valid:
            solution_count += 1
            
    # --- Output ---
    print("Problem: Find the number of ways to place non-attacking queens.")
    print(f"Board Size      = {BOARD_SIZE}x{BOARD_SIZE} (Toroidal)")
    print(f"Number of Queens = {NUM_QUEENS}")
    print("---------------------------------------")
    print(f"Result            = {solution_count}")

# Execute the function to find and print the solution
solve_toroidal_queens()