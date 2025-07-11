import itertools

def count_non_attacking_queens():
    """
    Calculates the number of ways to place a specified number of non-attacking 
    queens on a toroidal chessboard of a given size.
    """
    # The size of the toroidal board
    board_size = 5
    # The number of non-attacking queens to place
    num_queens = 4

    def are_queens_attacking(q1, q2):
        """
        Checks if two queens on a toroidal board are attacking each other.
        A board of size N x N is assumed from the board_size variable.
        q1 and q2 are tuples representing the (row, col) coordinates of the queens.
        """
        r1, c1 = q1
        r2, c2 = q2
        
        # Check for same row or same column
        if r1 == r2 or c1 == c2:
            return True
            
        # Check for same main diagonal (toroidal)
        if (r1 - c1) % board_size == (r2 - c2) % board_size:
            return True
            
        # Check for same anti-diagonal (toroidal)
        if (r1 + c1) % board_size == (r2 + c2) % board_size:
            return True
            
        return False

    # Generate a list of all squares on the board
    all_squares = [(r, c) for r in range(board_size) for c in range(board_size)]
    
    # Generate all possible combinations for placing the queens
    all_placements = itertools.combinations(all_squares, num_queens)
    
    valid_placements_count = 0
    # Iterate through each placement combination
    for placement in all_placements:
        is_valid = True
        # Check all pairs of queens within the current placement
        for queen1, queen2 in itertools.combinations(placement, 2):
            if are_queens_attacking(queen1, queen2):
                is_valid = False
                break  # Move to the next placement
        
        if is_valid:
            valid_placements_count += 1
            
    # The problem is about placing 4 queens on a 5x5 board.
    # The final equation is simply the total count.
    print(f"The number of ways to place {num_queens} non-attacking queens on a {board_size}x{board_size} toroidal chessboard is:")
    print(valid_placements_count)

if __name__ == '__main__':
    count_non_attacking_queens()