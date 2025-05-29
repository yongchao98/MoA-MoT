def print_board(board):
    print("   A B C D E F G H I J K L M")
    for i in range(13, 0, -1):
        print(f"{i:2d}", end=" ")
        for j in range(13):
            print(board[i-1][j], end=" ")
        print()

def count_liberties(board, x, y, color, visited=None):
    if visited is None:
        visited = set()
    
    if x < 0 or x >= 13 or y < 0 or y >= 13:
        return 0
    
    if board[y][x] != color:
        if board[y][x] == '.':
            return 1
        return 0
    
    if (x, y) in visited:
        return 0
        
    visited.add((x, y))
    liberties = 0
    
    # Check all adjacent positions
    for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
        new_x, new_y = x + dx, y + dy
        liberties += count_liberties(board, new_x, new_y, color, visited)
    
    return liberties

# Initialize the board
board = [['.'] * 13 for _ in range(13)]

# Set up the relevant stones
stones = {
    'X': [(2, 12), (2, 9), (8, 9), (2, 4), (2, 3), (1, 2), (0, 1), (3, 0), (2, 3), (3, 2), (10, 6), (10, 10)],
    'O': [(7, 12), (0, 9), (2, 2), (2, 1), (2, 0), (4, 2), (12, 8)]
}

for color in stones:
    for x, y in stones[color]:
        board[y][x] = color

print_board(board)

# Check liberties of white group in corner
white_group_liberties = count_liberties(board, 2, 2, 'O')
print(f"\nLiberties of white group: {white_group_liberties}")

# Check potential vital points
vital_points = [(1, 0), (1, 1), (2, 1)]
for x, y in vital_points:
    if board[y][x] == '.':
        print(f"Checking point ({chr(65+x)}{y+1})")