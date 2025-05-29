def print_board(board):
    for y in range(10, -1, -1):
        row = []
        for x in range(11):
            row.append(board[y][x])
        print(f"{y:2d} {''.join(row)}")
    print("   ABCDEFGHIJK")

# Create the board (11x11)
board = [
    ['.'] * 11 for _ in range(11)
]

# Place the stones as given in the problem
# Black stones (X)
black_positions = [
    (0,10), (5,10), (6,10), (7,10),  # Row 11
    (8,9),  # Row 10
    (3,8),  # Row 9
    (4,7), (6,7),  # Row 8
    (5,4),  # Row 5
    (9,3),  # Row 4
    (9,2),  # Row 3
]

# White stones (O)
white_positions = [
    (6,9), (7,9),  # Row 10
    (6,8),  # Row 9
    (3,4),  # Row 5
    (8,4),  # Row 5
    (2,1),  # Row 2
]

# Place stones on the board
for x, y in black_positions:
    board[y][x] = 'X'
for x, y in white_positions:
    board[y][x] = 'O'

# Check if H9 (7,8) would capture
def check_capture(board, x, y):
    # Place a black stone at the target position
    board[y][x] = 'X'
    
    # Check if any adjacent white stones are captured
    captured = False
    adjacent_positions = [(7,9), (6,8), (6,9)]  # Positions of the white stones we're trying to capture
    
    # Check if these white stones have any liberties left
    for wx, wy in adjacent_positions:
        has_liberty = False
        # Check all adjacent positions
        for dx, dy in [(0,1), (0,-1), (1,0), (-1,0)]:
            nx, ny = wx + dx, wy + dy
            if 0 <= nx < 11 and 0 <= ny < 11:
                if board[ny][nx] == '.':
                    has_liberty = True
                    break
        if not has_liberty:
            captured = True
            break
    
    # Reset the board
    board[y][x] = '.'
    return captured

# Test H9 move (7,8)
result = check_capture(board, 7, 8)
print(f"Can capture at H9: {result}")
print("\nBoard position:")
print_board(board)