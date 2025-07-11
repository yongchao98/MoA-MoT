import itertools

def solve_toroidal_queens_problem():
    """
    Calculates the number of ways to place a specified number of non-attacking 
    queens on a toroidal chessboard of a given size.
    """
    board_dim = 5
    num_queens = 4

    def is_attacking(p1, p2, n):
        """
        Checks if two queens at positions p1 and p2 on an n x n toroidal board 
        are attacking each other.
        """
        r1, c1 = p1
        r2, c2 = p2
        # Check for row or column attack
        if r1 == r2 or c1 == c2:
            return True
        # Check for diagonal attacks using modular arithmetic for the wrap-around
        if (r1 - c1) % n == (r2 - c2) % n:
            return True
        if (r1 + c1) % n == (r2 + c2) % n:
            return True
        return False

    # A list of all coordinates on the board
    squares = [(r, c) for r in range(board_dim) for c in range(board_dim)]

    # Generate all combinations of squares to place the queens
    possible_placements = itertools.combinations(squares, num_queens)

    valid_count = 0
    # Iterate through each possible placement
    for placement in possible_placements:
        is_valid_placement = True
        # Check all pairs of queens in the current placement
        for queen_pair in itertools.combinations(placement, 2):
            if is_attacking(queen_pair[0], queen_pair[1], board_dim):
                is_valid_placement = False
                break  # This placement is invalid, so we stop checking its pairs
        
        if is_valid_placement:
            valid_count += 1
            
    # The final "equation" is the summary of the problem and its solution,
    # containing all relevant numbers.
    print(f"On a {board_dim}x{board_dim} toroidal chessboard, the number of ways to place {num_queens} non-attacking queens is {valid_count}.")

solve_toroidal_queens_problem()