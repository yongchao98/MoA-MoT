# Define the board
board = [
    ['.', '.', '.', '.', 'X', '.', 'X', '.', '.', '.', '.'],
    ['.', 'X', '.', 'X', '.', '.', '.', 'X', 'X', '.', '.'],
    ['X', 'O', 'O', 'O', 'X', '.', '.', 'X', '.', '.', '.'],
    ['X', 'X', 'O', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', 'O', '.'],
    ['.', '.', '.', '.', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'O', '.', '.', '.', '.', 'O', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Function to count liberties of a group
def count_liberties(board, x, y):
    visited = set()
    stack = [(x, y)]
    liberties = 0
    color = board[x][y]
    
    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))
        
        # Check all 4 directions
        for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nx, ny = cx + dx, cy + dy
            if 0 <= nx < len(board) and 0 <= ny < len(board[0]):
                if board[nx][ny] == '.':
                    liberties += 1
                elif board[nx][ny] == color and (nx, ny) not in visited:
                    stack.append((nx, ny))
    
    return liberties

# Simulate the move D8
board[3][3] = 'X'

# Count liberties of the White group at C9, D9, E9
liberties_C9 = count_liberties(board, 2, 1)
liberties_D9 = count_liberties(board, 2, 2)
liberties_E9 = count_liberties(board, 2, 3)

print(liberties_C9, liberties_D9, liberties_E9)