def is_safe(queens, row, col):
    for queen_row, queen_col in queens:
        if (queen_row == row or  # Same row
            queen_col == col or  # Same column
            abs(queen_row - row) == abs(queen_col - col)):  # Same diagonal
            return False
    return True

def solve_queens():
    # Start with fixed queens
    queens = [(0,1), (3,6), (5,7)]
    
    def backtrack(row):
        if len(queens) == 8:
            return True
            
        # Skip rows that already have queens
        while row < 8 and any(q[0] == row for q in queens):
            row += 1
            
        if row >= 8:
            return False
            
        for col in range(8):
            if is_safe(queens, row, col):
                queens.append((row, col))
                if backtrack(row + 1):
                    return True
                queens.pop()
        return False
    
    if backtrack(0):
        # Sort positions for consistent output
        queens.sort()
        result = ", ".join(f"{pos[0]} {pos[1]}" for pos in queens)
        print(f"<<<{result}>>>")
    else:
        print("No solution exists")

solve_queens()