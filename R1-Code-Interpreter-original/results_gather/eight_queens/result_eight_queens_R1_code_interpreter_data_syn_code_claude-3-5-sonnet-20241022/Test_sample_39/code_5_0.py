def is_safe(queens, row, col, forbidden):
    if (row, col) in forbidden:
        return False
    
    # Check against each existing queen
    for qr, qc in queens:
        # Same row, column, or diagonal check
        if (row == qr) or (col == qc) or (abs(row - qr) == abs(col - qc)):
            return False
    return True

def solve_queens_puzzle():
    # Initialize with pre-placed queen and forbidden positions
    queens = [(0, 5)]  # Pre-placed queen
    forbidden = {(0, 6), (3, 0), (5, 0), (7, 3)}
    
    def solve_column(col, remaining):
        if remaining == 0:
            return True
            
        # Skip column 5 as it already has a queen
        if col == 5:
            return solve_column(col + 1, remaining)
            
        if col >= 8:
            return False
            
        # Try each row in current column
        for row in range(8):
            if is_safe(queens, row, col, forbidden):
                queens.append((row, col))
                if solve_column(col + 1, remaining - 1):
                    return True
                queens.pop()
                
        return False
    
    # Start solving from column 0, need to place 7 more queens
    if solve_column(0, 7):
        # Sort queens by column for consistent output
        queens.sort(key=lambda x: x[1])
        print(",".join(f"{row} {col}" for row, col in queens))
        return True
    return False

solve_queens_puzzle()