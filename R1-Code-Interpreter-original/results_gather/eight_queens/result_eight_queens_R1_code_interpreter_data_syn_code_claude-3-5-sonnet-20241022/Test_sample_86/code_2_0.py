def is_attacking(pos1, pos2):
    # Check if two positions are attacking each other
    row1, col1 = pos1
    row2, col2 = pos2
    return (row1 == row2 or  # same row
            col1 == col2 or  # same column
            abs(row1 - row2) == abs(col1 - col2))  # same diagonal

def is_valid_position(position, queens, blocked_pos):
    row, col = position
    # Check if position is blocked
    if position == blocked_pos:
        return False
    # Check if position conflicts with any existing queen
    return not any(is_attacking(position, queen) for queen in queens)

def solve_queens(n, fixed_queen, blocked_pos):
    queens = {fixed_queen}  # Start with the fixed queen
    rows = set(range(n))
    cols = set(range(n))
    
    def backtrack(queens):
        if len(queens) == n:
            return queens
        
        # Get all remaining positions to try
        for row in range(n):
            for col in range(n):
                if (row, col) in queens:
                    continue
                    
                pos = (row, col)
                if is_valid_position(pos, queens, blocked_pos):
                    queens.add(pos)
                    result = backtrack(queens)
                    if result:
                        return result
                    queens.remove(pos)
        return None

    result = backtrack(queens)
    return result

def find_queens_positions():
    n = 8
    fixed_queen = (2, 6)  # Pre-placed queen
    blocked_pos = (6, 3)  # Blocked position
    
    solution = solve_queens(n, fixed_queen, blocked_pos)
    if solution:
        # Convert to required format and sort
        positions = [f"{row} {col}" for row, col in sorted(solution)]
        print("<<<" + ", ".join(positions) + ">>>")
    else:
        print("No solution exists")

find_queens_positions()