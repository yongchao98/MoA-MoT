def print_board(board):
    print("   A B C D E F G H I J K L")
    for i in range(12, 0, -1):
        print(f"{i:2d}", end=" ")
        for j in range(12):
            print(board[i-1][j], end=" ")
        print()

def get_liberties(board, x, y, visited):
    if x < 0 or x >= 12 or y < 0 or y >= 12:
        return set()
    if (x, y) in visited:
        return set()
    if board[y][x] == '.':
        return {(x, y)}
    if board[y][x] != 'O':
        return set()
    
    visited.add((x, y))
    liberties = set()
    for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
        liberties.update(get_liberties(board, x + dx, y + dy, visited))
    return liberties

# Initialize the board
board = [['.'] * 12 for _ in range(12)]
# Fill in the stones from the given position
stones = """
. . . . . O . . . . . O
. . . . . O . . . . . .
. . . . . . . X . . O X
. . . . . O X . X . . .
. O . . . X O O O X . O
X . . . . . X O X . . .
. . . . . . . X . . . .
. . . . . . . . . . . .
. . . . . . O . . . . .
. . . . . . . . . . . .
. . . . . . . . . . . .
. . . . . . . . . . . .
"""

# Parse the stones string into the board
for i, line in enumerate(stones.strip().split('\n')):
    for j, char in enumerate(line.split()):
        board[11-i][j] = char

# Check liberties of the central white group
central_group_liberties = get_liberties(board, 7, 6, set())  # H7 position
print("Liberties of central white group:", end=" ")
for x, y in central_group_liberties:
    print(f"({chr(65+x)}{y+1})", end=" ")
print()