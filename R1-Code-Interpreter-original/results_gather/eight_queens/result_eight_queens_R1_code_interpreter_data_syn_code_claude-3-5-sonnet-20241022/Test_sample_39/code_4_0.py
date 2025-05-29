def is_valid(queens, row, col, forbidden):
    # Check if position is forbidden
    if (row, col) in forbidden:
        return False
        
    # Check against each existing queen
    for queen_row, queen_col in queens:
        if queen_row == row:  # Same row
            return False
        if queen_col == col:  # Same column
            return False
        if abs(queen_row - row) == abs(queen_col - col):  # Diagonal
            return False
    return True

def solve_queens():
    # Initialize constraints
    forbidden = {(0,6), (3,0), (5,0), (7,3)}
    # Start with the pre-placed queen
    queens = [(0,5)]
    
    def backtrack(row):
        if len(queens) == 8:
            return True
            
        # Move to next row if current row has a queen
        while row < 8 and any(q[0] == row for q in queens):
            row += 1
            
        if row >= 8:
            return False
            
        # Try each column in current row
        for col in range(8):
            if is_valid(queens, row, col, forbidden):
                queens.append((row, col))
                if backtrack(row + 1):
                    return True
                queens.pop()
        return False
    
    # Start with row 1 since row 0 has pre-placed queen
    if backtrack(1):
        # Sort queens by row and format output
        queens.sort()
        print(",".join(f"{row} {col}" for row, col in queens))
    else:
        print("No solution exists")

solve_queens()