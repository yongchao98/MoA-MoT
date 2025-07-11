# The size of the board
N = 5
# The number of queens to place
Q = 4

# A global variable to store the total count of solutions
solution_count = 0

def is_safe(r, c, placements):
    """
    Checks if placing a queen at (r, c) is safe with respect to the queens
    already in 'placements'.
    On a toroidal board, attacks are checked using modulo arithmetic.
    """
    for pr, pc in placements:
        # Check for same row or same column (this check is sufficient for standard boards)
        if r == pr or c == pc:
            return False
        # Check for same diagonal, which wraps around the board
        if (r - c) % N == (pr - pc) % N:
            return False
        # Check for same anti-diagonal, which also wraps around
        if (r + c) % N == (pr + pc) % N:
            return False
    return True

def solve_toroidal_queens(k, start_index, placements):
    """
    A recursive function to find solutions using backtracking.
    k: The number of queens already placed.
    start_index: The linear index on the board to start searching from. This ensures
                 that we find unique combinations of squares, not permutations.
    placements: A list of (row, col) tuples for the queens placed so far.
    """
    global solution_count

    # If all queens have been placed, we have found a valid solution.
    if k == Q:
        solution_count += 1
        return

    # Iterate through all squares on the board, starting from 'start_index'.
    for i in range(start_index, N * N):
        r = i // N
        c = i % N

        # If it's safe to place a queen at the current square (r, c)...
        if is_safe(r, c, placements):
            # ...place the queen and recurse to place the next one.
            placements.append((r, c))
            solve_toroidal_queens(k + 1, i + 1, placements)
            # Backtrack: remove the queen to explore other possibilities.
            placements.pop()

def main():
    """
    Main function to initiate the solver and print the result.
    """
    # Start the backtracking process.
    # We begin with 0 queens placed, starting the search from square 0.
    solve_toroidal_queens(0, 0, [])

    # Print the final count of solutions.
    print(f"The number of ways to place {Q} non-attacking queens on a {N}x{N} toroidal chessboard is: {solution_count}")

if __name__ == "__main__":
    main()
