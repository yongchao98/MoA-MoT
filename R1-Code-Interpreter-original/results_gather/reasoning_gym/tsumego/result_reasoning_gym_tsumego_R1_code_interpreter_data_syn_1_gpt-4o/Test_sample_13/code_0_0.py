# Define the board
board = [
    ['.', '.', '.', '.', '.', 'O', 'X', 'X', '.', '.', '.'],
    ['.', '.', 'X', '.', 'O', 'X', 'O', 'O', 'X', '.', '.'],
    ['O', '.', '.', 'O', '.', 'X', 'O', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', '.', '.', '.', '.', '.'],
    ['.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'O', '.', '.', '.', '.', '.', 'X', '.', '.', '.'],
    ['X', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['O', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Function to check liberties of a group
def check_liberties(board, x, y, color):
    visited = set()
    stack = [(x, y)]
    liberties = 0
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    
    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))
        
        for dx, dy in directions:
            nx, ny = cx + dx, cy + dy
            if 0 <= nx < len(board) and 0 <= ny < len(board[0]):
                if board[nx][ny] == '.':
                    liberties += 1
                elif board[nx][ny] == color and (nx, ny) not in visited:
                    stack.append((nx, ny))
    
    return liberties

# Check liberties for the group at (E10, F10, G10, F9)
liberties = check_liberties(board, 1, 4, 'O')
print("Liberties for the group at (E10, F10, G10, F9):", liberties)

# Suggest move
suggested_move = 'H10'
print("Suggested move for Black:", suggested_move)