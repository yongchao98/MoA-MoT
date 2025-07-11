import itertools

def solve_toroidal_queens():
    """
    This function calculates the number of ways to place a specified number of 
    non-attacking queens on a toroidal chessboard.
    """
    n = 5  # The size of the board's side (e.g., 5 for a 5x5 board)
    k = 4  # The number of queens to place

    def is_attacking(q1, q2):
        """
        Checks if two queens, q1 and q2, on an nxn toroidal board are attacking each other.
        Each queen is represented by a tuple (row, col).
        """
        r1, c1 = q1
        r2, c2 = q2

        # Check for attack in the same row or column
        if r1 == r2 or c1 == c2:
            return True
        
        # Check for attack on the toroidal diagonals
        # Main diagonal condition: (r - c) mod n is constant
        if (r1 - c1) % n == (r2 - c2) % n:
            return True
            
        # Anti-diagonal condition: (r + c) mod n is constant
        if (r1 + c1) % n == (r2 + c2) % n:
            return True
            
        return False

    # Generate all possible squares on the board as (row, col) tuples
    all_squares = [(r, c) for r in range(n) for c in range(n)]
    
    # Initialize a counter for the number of valid placements found
    valid_placements_count = 0
    
    # Iterate through all unique combinations of k squares to place the queens
    for placement in itertools.combinations(all_squares, k):
        is_valid = True
        # For each placement, check every pair of queens to see if they attack each other
        for queen1, queen2 in itertools.combinations(placement, 2):
            if is_attacking(queen1, queen2):
                is_valid = False
                break  # If a single pair attacks, the placement is invalid
        
        if is_valid:
            valid_placements_count += 1
            
    # Print the final result in a descriptive sentence
    print(f"On a {n}x{n} toroidal chessboard, the number of ways to place {k} non-attacking queens is {valid_placements_count}.")

# Execute the function to find and print the solution
solve_toroidal_queens()