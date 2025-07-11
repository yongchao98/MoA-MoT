import sys

def solve():
    """
    Solves the beautiful matrix problem.
    """
    try:
        # Read the matrix size N
        n = int(sys.stdin.readline())
    except (ValueError, IndexError):
        # Handle cases with no or invalid input for N
        return

    k_val = 0
    k_row = 0
    k_col = 0

    # Read the matrix to find the coordinates of the non-zero digit
    # Rows are 1-indexed from 1 to N
    for i in range(1, n + 1):
        try:
            line = sys.stdin.readline().split()
            # Columns are 1-indexed from 1 to N
            for j in range(1, n + 1):
                val = int(line[j - 1])
                if val != 0:
                    k_val = val
                    k_row = i
                    k_col = j
                    # Since there's only one non-zero number, we could break
                    # here for efficiency, but completing the loops is fine.
        except (ValueError, IndexError):
            # Skip malformed lines
            continue
    
    # The center of the matrix (1-indexed)
    center = n // 2 + 1

    # The number of moves is the Manhattan distance from the current position to the center
    # This is the "final equation" for the number of moves, z.
    # z = |k_row - center| + |k_col - center|
    moves = abs(k_row - center) + abs(k_col - center)
    
    # Print the integers k, r, c, and z as required.
    # We are outputting k, r, c which are inputs to the equation, and z which is the result.
    print(k_val, k_row, k_col, moves)

solve()