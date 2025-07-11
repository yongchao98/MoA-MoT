# Set the size of the board and the number of queens.
N = 5
NUM_QUEENS = 4

# This global variable will store the total count of valid solutions.
solution_count = 0

def is_safe(new_pos, queens):
    """
    Checks if a new queen at new_pos is not attacked by any queen in the queens list.
    The check is performed on a toroidal (wrap-around) board.
    
    Args:
        new_pos: A tuple (row, col) for the new queen.
        queens: A list of (row, col) tuples for existing queens.

    Returns:
        True if the position is safe, False otherwise.
    """
    r2, c2 = new_pos
    for r1, c1 in queens:
        # Check for same row or column (standard attack).
        if r1 == r2 or c1 == c2:
            return False
        # Check for toroidal diagonals.
        # Two queens are on the same diagonal if (r1-c1) % N == (r2-c2) % N.
        if (r1 - c1) % N == (r2 - c2) % N:
            return False
        # Check for toroidal anti-diagonals.
        # Two queens are on the same anti-diagonal if (r1+c1) % N == (r2+c2) % N.
        if (r1 + c1) % N == (r2 + c2) % N:
            return False
    return True

def solve_placements(queens_to_place, start_index, placed_queens):
    """
    Uses backtracking to find all valid placements for the queens.
    
    Args:
        queens_to_place: The number of queens remaining to be placed.
        start_index: The linear board index to start searching from. This ensures
                     we count unique combinations of squares, not permutations.
        placed_queens: A list of coordinates for queens already placed.
    """
    global solution_count

    # Base case: If no more queens need to be placed, we have found a valid solution.
    if queens_to_place == 0:
        solution_count += 1
        return

    # Iterate through all board squares, starting from the given index.
    for i in range(start_index, N * N):
        row = i // N
        col = i % N
        
        # Check if placing a queen at the current (row, col) is safe.
        if is_safe((row, col), placed_queens):
            # If safe, place the queen (add it to our list).
            placed_queens.append((row, col))
            # Recursively call the function to place the next queen.
            # The next search starts from the next index to avoid duplicates.
            solve_placements(queens_to_place - 1, i + 1, placed_queens)
            # Backtrack: remove the queen to explore other possibilities.
            placed_queens.pop()

# Initial call to start the search process.
# We need to place NUM_QUEENS, starting the search from index 0, with an empty board.
solve_placements(NUM_QUEENS, 0, [])

# The final "equation" showing all the numbers involved in the problem.
print(f"Number of ways to place {NUM_QUEENS} non-attacking queens on a {N}x{N} toroidal board = {solution_count}")