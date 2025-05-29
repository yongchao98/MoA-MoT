def count_liberties(board, x, y):
    """Count the liberties of a group of stones starting from (x, y)."""
    visited = set()
    stack = [(x, y)]
    liberties = 0
    color = board[x][y]
    
    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))
        
        # Check all four directions
        for nx, ny in [(cx-1, cy), (cx+1, cy), (cx, cy-1), (cx, cy+1)]:
            if 0 <= nx < len(board) and 0 <= ny < len(board[0]):
                if board[nx][ny] == '.':
                    liberties += 1
                elif board[nx][ny] == color and (nx, ny) not in visited:
                    stack.append((nx, ny))
    
    return liberties

# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'O', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', 'O', 'X', '.', '.', '.', '.', '.', 'O', '.', '.'],
    ['X', 'O', 'O', '.', '.', '.', '.', 'O', '.', '.', '.', '.'],
    ['.', 'X', 'O', 'X', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', '.', 'X', 'X'],
    ['.', '.', '.', '.', 'O', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['X', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', '.', 'X', 'O', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Evaluate potential moves
liberties_C9 = count_liberties(board, 3, 2)  # D9
liberties_D8 = count_liberties(board, 4, 2)  # C8

print("Liberties after playing at C9:", liberties_C9)
print("Liberties after playing at D8:", liberties_D8)