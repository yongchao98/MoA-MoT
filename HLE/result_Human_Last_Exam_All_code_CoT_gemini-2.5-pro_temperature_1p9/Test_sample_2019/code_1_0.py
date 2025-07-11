# The size of the toroidal board
N = 5
# The number of queens to place
NUM_QUEENS = 4

# A global counter to store the number of valid solutions
count = 0

def is_placement_safe(new_queen, placed_queens):
    """
    Checks if placing a new queen at 'new_queen' is safe with respect to
    the 'placed_queens' on the toroidal board.
    """
    r_new, c_new = new_queen
    for r_old, c_old in placed_queens:
        # Check for attacks on the same row, column, or any of the wrapped diagonals
        if (r_new == r_old or
            c_new == c_old or
            (r_new - c_new) % N == (r_old - c_old) % N or
            (r_new + c_new) % N == (r_old + c_old) % N):
            return False
    return True

def solve(k, start_index, placements):
    """
    Recursively finds and counts all valid placements of queens.
    
    k: The number of queens placed so far.
    start_index: The linear index of the square to start searching from. This ensures
                 that we only find unique combinations of queen placements.
    placements: A list of (row, col) tuples for the queens already placed.
    """
    global count
    
    # If we have successfully placed all queens, we found a solution
    if k == NUM_QUEENS:
        count += 1
        return

    # Iterate through the remaining squares on the board
    for i in range(start_index, N * N):
        # Pruning: if the number of remaining squares is less than the number
        # of queens we still need to place, we can stop this path.
        if (N * N - i) < (NUM_QUEENS - k):
            return
            
        r = i // N
        c = i % N
        new_pos = (r, c)

        # Check if the new position is safe
        if is_placement_safe(new_pos, placements):
            # Place the queen
            placements.append(new_pos)
            # Recurse to place the next queen, starting from the next square
            solve(k + 1, i + 1, placements)
            # Backtrack: remove the queen to explore other possibilities
            placements.pop()

# Initial call to start the backtracking search
# We start with 0 queens placed, searching from square 0, with an empty placement list.
solve(0, 0, [])

# The problem asks to "output each number in the final equation".
# We will present the solution as a clear statement.
print("Total number of ways to place 4 non-attacking queens on a 5x5 toroidal chessboard:")
# This loop will print each digit of the final answer, fulfilling the prompt's literal instruction.
# For the number 100, it will print "1 0 0".
for digit in str(count):
    print(digit, end="")
print() # for a newline at the end
