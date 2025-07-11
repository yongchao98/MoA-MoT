import sys

# Suppress recursion limit errors for deep searches, though not strictly
# necessary for this small problem size.
sys.setrecursionlimit(2000)

# Board configuration
N = 5  # Board size is NxN
K = 4  # Number of queens to place

def is_safe(r_new, c_new, placed_queens):
    """
    Checks if a new queen at (r_new, c_new) can be placed without
    being attacked by any of the queens in placed_queens.
    This check includes toroidal (wrapped) diagonal attacks.
    """
    for r_old, c_old in placed_queens:
        # Check for same row or column (should not happen with the way
        # we construct solutions but included for clarity).
        if r_new == r_old or c_new == c_old:
            return False
            
        # Check for toroidal main diagonals.
        # This condition checks if the queens lie on the same "r-c" diagonal.
        if (r_new - c_new + N) % N == (r_old - c_old + N) % N:
            return False
            
        # Check for toroidal anti-diagonals.
        # This condition checks if the queens lie on the same "r+c" diagonal.
        if (r_new + c_new) % N == (r_old + c_old) % N:
            return False
            
    return True

def count_placements(queens_to_place, start_index, placed_queens):
    """
    Recursively finds the number of ways to place the remaining queens.
    
    :param queens_to_place: The number of queens we still need to place.
    :param start_index: The linear cell index to start searching from. This
                        prevents counting permutations of the same placement.
    :param placed_queens: A list of (row, col) tuples for queens already on the board.
    :return: The total count of valid placements from this state.
    """
    # Base case: If we have placed all K queens, we have found one valid solution.
    if queens_to_place == 0:
        return 1

    count = 0
    # Iterate through all possible cells starting from 'start_index'.
    for i in range(start_index, N * N):
        r = i // N
        c = i % N

        # Check if it's safe to place a queen here.
        if is_safe(r, c, placed_queens):
            # If safe, place the queen and recurse to place the next one.
            # We search for the next queen starting from the next cell (i + 1).
            new_placed_queens = placed_queens + [(r, c)]
            count += count_placements(queens_to_place - 1, i + 1, new_placed_queens)
    
    return count

# Start the backtracking search
# We need to place K queens, starting from cell 0, with an empty board.
total_ways = count_placements(K, 0, [])

# The final "equation" is the statement of the result.
# We print each number involved: N, K, and the final count.
print(f"On a {N}x{N} toroidal board, {K} non-attacking queens can be placed in {total_ways} ways.")

<<<250>>>