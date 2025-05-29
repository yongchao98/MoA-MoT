def is_valid(positions, row, col):
    # Check if there's a queen in the same column or diagonals
    for i in range(row):
        if positions[i] == col or \
           positions[i] - i == col - row or \
           positions[i] + i == col + row:
            return False
    return True

def solve_n_queens(board, row, positions):
    if row == 8:
        return True
    
    for col in range(8):
        if board[row][col] == 'X':
            continue
        if is_valid(positions, row, col):
            positions[row] = col
            if solve_n_queens(board, row + 1, positions):
                return True
            positions[row] = -1
    
    return False

def place_queens():
    board = [
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 'X', 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 'X'],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 'X', 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0]
    ]
    
    positions = [-1] * 8
    if solve_n_queens(board, 0, positions):
        result = ', '.join(f"{i} {positions[i]}" for i in range(8))
        print(f"<<<{result}>>>")
    else:
        print("No solution found")

place_queens()