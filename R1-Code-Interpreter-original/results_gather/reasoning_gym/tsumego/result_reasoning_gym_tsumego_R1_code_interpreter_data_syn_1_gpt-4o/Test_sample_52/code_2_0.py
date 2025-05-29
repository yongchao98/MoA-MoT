# Define the board
board = [
    ['.', '.', 'X', '.', '.', '.', '.', 'X', '.', 'X', '.'],
    ['X', '.', '.', '.', 'X', '.', '.', '.', 'X', 'O', 'X'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', 'O', 'O', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', 'O', 'X'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', 'X', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'O', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'O', '.', 'X', '.', '.', '.']
]

# Function to check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 11 and 0 <= y < 11

# Function to get the liberties of a group of stones
def get_liberties(x, y, color):
    visited = set()
    stack = [(x, y)]
    liberties = set()
    
    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))
        
        # Check all adjacent positions
        for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nx, ny = cx + dx, cy + dy
            if is_within_board(nx, ny):
                if board[nx][ny] == '.':
                    liberties.add((nx, ny))
                elif board[nx][ny] == color and (nx, ny) not in visited:
                    stack.append((nx, ny))
    
    return liberties

# Check the liberties of the white group at J10, I9, J9
white_group_liberties = get_liberties(1, 9, 'O')

# Check if playing at J8 captures the group
if (3, 9) in white_group_liberties:
    white_group_liberties.remove((3, 9))

# If no liberties remain, the group is captured
captured = len(white_group_liberties) == 0

print("Captured:", captured)