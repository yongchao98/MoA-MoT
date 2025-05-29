# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'O', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'O', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['X', 'O', '.', '.', 'O', '.', '.', '.', '.', '.', '.', '.'],
    ['X', 'O', 'O', 'X', '.', '.', '.', 'X', '.', '.', 'X', '.'],
    ['O', 'X', 'X', '.', '.', '.', '.', '.', 'O', '.', 'O', 'X'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Function to check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 12 and 0 <= y < 12

# Function to get the liberties of a group
def get_liberties(x, y, color):
    visited = set()
    stack = [(x, y)]
    liberties = set()
    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))
        for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nx, ny = cx + dx, cy + dy
            if is_within_board(nx, ny):
                if board[nx][ny] == '.':
                    liberties.add((nx, ny))
                elif board[nx][ny] == color and (nx, ny) not in visited:
                    stack.append((nx, ny))
    return liberties

# Check the group at (10, 9) and (10, 11)
group_liberties = get_liberties(10, 9, 'O')
if len(group_liberties) == 1:
    move = group_liberties.pop()
    print(f"<<<{chr(move[1] + ord('A'))}{12 - move[0]}>>>")
else:
    print("No single move capture available.")