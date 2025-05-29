def get_neighbors(x, y, size=9):
    neighbors = []
    for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < size and 0 <= new_y < size:
            neighbors.append((new_x, new_y))
    return neighbors

def find_group(board, x, y, color, visited=None):
    if visited is None:
        visited = set()
    
    if (x, y) in visited or board[y][x] != color:
        return visited
    
    visited.add((x, y))
    for nx, ny in get_neighbors(x, y):
        find_group(board, nx, ny, color, visited)
    return visited

def count_liberties(board, group):
    liberties = set()
    for x, y in group:
        for nx, ny in get_neighbors(x, y):
            if board[ny][nx] == '.':
                liberties.add((nx, ny))
    return liberties

# Initialize the board
board = [
    ['.', 'X', '.', '.', '.', '.', '.', '.', '.'],  # 9
    ['O', '.', 'X', '.', '.', '.', '.', '.', '.'],  # 8
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.'],  # 7
    ['.', '.', 'X', '.', 'X', '.', '.', 'X', '.'],  # 6
    ['.', 'X', 'O', 'O', 'O', 'X', '.', '.', '.'],  # 5
    ['O', '.', 'X', 'O', 'X', '.', '.', '.', '.'],  # 4
    ['O', '.', '.', 'X', '.', '.', '.', '.', '.'],  # 3
    ['.', '.', '.', 'O', '.', 'O', '.', '.', '.'],  # 2
    ['.', '.', 'O', '.', '.', 'X', '.', '.', '.']   # 1
]

# Find all white groups and their liberties
white_groups = []
visited_whites = set()

for y in range(9):
    for x in range(9):
        if board[y][x] == 'O' and (x, y) not in visited_whites:
            group = find_group(board, x, y, 'O')
            if group:
                visited_whites.update(group)
                liberties = count_liberties(board, group)
                white_groups.append((group, liberties))

# Analyze each white group
print("White groups analysis:")
for i, (group, liberties) in enumerate(white_groups):
    print(f"Group {i+1}: Stones at {sorted(group)}")
    print(f"Liberties at {sorted(liberties)}")
    print(f"Number of liberties: {len(liberties)}\n")