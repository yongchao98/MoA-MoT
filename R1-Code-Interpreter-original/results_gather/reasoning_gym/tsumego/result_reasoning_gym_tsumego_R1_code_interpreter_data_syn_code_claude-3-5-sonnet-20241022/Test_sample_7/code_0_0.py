def is_valid(x, y):
    return 0 <= x < 10 and 0 <= y < 10

def get_liberties(board, x, y, visited):
    if not is_valid(x, y) or (x, y) in visited:
        return set()
    
    if board[y][x] == '.':
        return {(x, y)}
    
    if board[y][x] != board[y0][x0]:
        return set()
        
    visited.add((x, y))
    liberties = set()
    for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
        new_x, new_y = x + dx, y + dy
        liberties.update(get_liberties(board, new_x, new_y, visited))
    return liberties

# Create board representation
board = [
    ['.', '.', '.', 'O', 'O', '.', '.', '.', '.', '.'],  # 10
    ['.', '.', '.', '.', '.', 'X', 'O', '.', '.', '.'],  # 9
    ['X', '.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 8
    ['.', '.', '.', '.', '.', '.', 'O', '.', '.', '.'],  # 7
    ['.', '.', '.', '.', 'O', '.', 'X', '.', '.', '.'],  # 6
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X'],  # 5
    ['.', '.', '.', '.', '.', 'X', '.', '.', '.', '.'],  # 4
    ['.', '.', '.', '.', '.', '.', 'X', 'X', 'X', '.'],  # 3
    ['.', '.', '.', '.', '.', 'X', 'X', 'O', 'O', 'X'],  # 2
    ['.', '.', '.', '.', '.', '.', 'X', 'O', '.', '.']   # 1
]

# Check the bottom right group
x0, y0 = 7, 9  # H1
visited = set()
liberties = get_liberties(board, x0, y0, visited)
print(f"Liberties for white group at H1: {liberties}")

# The key liberty point would be at I1 (8, 9)
if (8, 9) in liberties:
    print("Playing at I1 would capture the white group")