import itertools

def count_non_attacking_queens():
    """
    Calculates the number of ways to place 4 non-attacking queens on a 5x5 toroidal chessboard.
    """
    N = 5  # Size of the board
    K = 4  # Number of queens

    def is_attacking(q1, q2):
        """Checks if two queens q1 and q2 attack each other on an N x N toroidal board."""
        r1, c1 = q1
        r2, c2 = q2

        # Check for attacks on the same row, column, or diagonals
        if r1 == r2 or c1 == c2:
            return True
        if (r1 - c1) % N == (r2 - c2) % N:
            return True
        if (r1 + c1) % N == (r2 + c2) % N:
            return True
        return False

    squares = [(r, c) for r in range(N) for c in range(N)]
    valid_placements_count = 0

    # Iterate through all combinations of K squares to place queens
    for queens_pos in itertools.combinations(squares, K):
        is_valid = True
        # Check all pairs of queens in the current placement
        # There are C(K, 2) pairs to check.
        for i in range(K):
            for j in range(i + 1, K):
                if is_attacking(queens_pos[i], queens_pos[j]):
                    is_valid = False
                    break
            if not is_valid:
                break
        
        if is_valid:
            valid_placements_count += 1
    
    # The final equation is: Total Ways = valid_placements_count
    # The question asks to output each number in the final equation.
    # We will print the calculated number of ways.
    print(f"{valid_placements_count}")

count_non_attacking_queens()