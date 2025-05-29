def analyze_position(board, move_x, move_y):
    def is_inside(x, y):
        return 0 <= x < 13 and 0 <= y < 13
    
    def get_adjacent_points(x, y):
        points = []
        for dx, dy in [(0, 1), (0, -1), (1, 0), (-1, 0)]:
            new_x, new_y = x + dx, y + dy
            if is_inside(new_x, new_y):
                points.append((new_x, new_y))
        return points
    
    def count_liberties(x, y, color, visited=None):
        if visited is None:
            visited = set()
        if not is_inside(x, y) or board[y][x] != color:
            return 0
        pos = (x, y)
        if pos in visited:
            return 0
        visited.add(pos)
        liberties = 0
        for new_x, new_y in get_adjacent_points(x, y):
            if board[new_y][new_x] == '.':
                liberties += 1
            elif board[new_y][new_x] == color:
                liberties += count_liberties(new_x, new_y, color, visited)
        return liberties
    
    # Make a copy and play the move
    board_copy = [row[:] for row in board]
    board_copy[move_y][move_x] = 'X'
    
    # Check if the move creates any cuts
    cuts = 0
    vital_points = 0
    
    # Check surrounding points
    for x, y in get_adjacent_points(move_x, move_y):
        if is_inside(x, y) and board_copy[y][x] == 'O':
            liberties = count_liberties(x, y, 'O')
            if liberties <= 2:
                vital_points += 1
            if liberties == 1:
                cuts += 1
    
    return cuts, vital_points

# Create board
board = [['.'] * 13 for _ in range(13)]

# Set up stones
black_stones = [(1, 12), (3, 3), (2, 4), (0, 4), (2, 2), (5, 3)]
white_stones = [(1, 1), (1, 2), (1, 3), (2, 3), (0, 3)]

for x, y in black_stones:
    board[y][x] = 'X'
for x, y in white_stones:
    board[y][x] = 'O'

# Analyze key points
key_points = [
    (0, 1),  # A2
    (0, 2),  # A3
    (2, 1),  # C2
    (2, 0),  # C1
]

for x, y in key_points:
    cuts, vital = analyze_position(board, x, y)
    print(f"Move at ({chr(65+x)}{13-y}): Cuts={cuts}, Vital points={vital}")