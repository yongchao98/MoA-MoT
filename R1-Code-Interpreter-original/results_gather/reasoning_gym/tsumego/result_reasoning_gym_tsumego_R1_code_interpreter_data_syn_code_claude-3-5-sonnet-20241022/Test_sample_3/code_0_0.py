def get_liberties(board, x, y, color, visited=None):
    if visited is None:
        visited = set()
    
    if x < 0 or x >= 13 or y < 0 or y >= 13:
        return set()
    
    if board[y][x] != color:
        if board[y][x] == '.':
            return {(x, y)}
        return set()
    
    if (x, y) in visited:
        return set()
    
    visited.add((x, y))
    liberties = set()
    
    # Check all adjacent positions
    for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
        new_x, new_y = x + dx, y + dy
        liberties.update(get_liberties(board, new_x, new_y, color, visited))
    
    return liberties

# Create the board
board = [
    ['.', '.', '.', 'O', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'O', 'X', 'X', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', 'O', 'O', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', 'O', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.', 'X', '.', '.', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', '.', 'X', '.', 'X'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'O', '.', '.', '.', '.', '.', '.', '.'],
    ['X', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'O'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', '.', '.', 'O', '.', '.', '.', '.', '.', '.', 'O', '.']
]

# Check liberties of the white group at D11
liberties = get_liberties(board, 3, 2, 'O')  # D11 is at (3,2) in 0-based indexing
print(f"Liberties of white group at D11: {liberties}")