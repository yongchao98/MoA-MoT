from collections import deque

# Define the board and initial positions
board = {
    'A1': None, 'B1': 'w', 'C1': 'B', 'D1': None,
    'A2': None, 'B2': 'B', 'C2': None, 'D2': None,
    'A3': None, 'B3': None, 'C3': 'w', 'D3': None
}

# Define knight moves
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert board position to coordinates
def pos_to_coord(pos):
    col, row = pos
    return ord(col) - ord('A'), int(row) - 1

# Convert coordinates to board position
def coord_to_pos(x, y):
    return chr(x + ord('A')) + str(y + 1)

# Check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 4 and 0 <= y < 3

# Get possible moves for a knight
def get_knight_moves(pos):
    x, y = pos_to_coord(pos)
    moves = []
    for dx, dy in knight_moves:
        nx, ny = x + dx, y + dy
        if is_within_board(nx, ny):
            new_pos = coord_to_pos(nx, ny)
            if board[new_pos] is None:
                moves.append(new_pos)
    return moves

# BFS to find the sequence of moves
def bfs():
    queue = deque([('w', 'C3', 'B1', 'B2', 'C1', [])])
    visited = set()
    
    while queue:
        turn, w1, w2, b1, b2, path = queue.popleft()
        
        # Check if the goal is reached
        if (w1, w2) == ('B2', 'C1') and (b1, b2) == ('C3', 'B1'):
            return path
        
        # Generate next moves
        if turn == 'w':
            for knight in [w1, w2]:
                for move in get_knight_moves(knight):
                    new_w1, new_w2 = (move, w2) if knight == w1 else (w1, move)
                    if (new_w1, new_w2, b1, b2) not in visited:
                        visited.add((new_w1, new_w2, b1, b2))
                        queue.append(('B', new_w1, new_w2, b1, b2, path + [f"w,{knight},{move}"]))
        else:
            for knight in [b1, b2]:
                for move in get_knight_moves(knight):
                    new_b1, new_b2 = (move, b2) if knight == b1 else (b1, move)
                    if (w1, w2, new_b1, new_b2) not in visited:
                        visited.add((w1, w2, new_b1, new_b2))
                        queue.append(('w', w1, w2, new_b1, new_b2, path + [f"B,{knight},{move}"]))
    
    return "No"

# Find the solution
solution = bfs()
print(solution)