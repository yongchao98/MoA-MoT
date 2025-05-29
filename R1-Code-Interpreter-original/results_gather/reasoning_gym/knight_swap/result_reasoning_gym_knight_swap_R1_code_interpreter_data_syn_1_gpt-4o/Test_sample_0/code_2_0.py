from collections import deque

# Define the board and initial positions
board = {
    'A1': None, 'B1': None, 'C1': 'w', 'D1': None,
    'A2': 'B', 'B2': None, 'C2': None, 'D2': None,
    'A3': None, 'B3': 'w', 'C3': 'B', 'D3': None
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

# Find all possible moves for a knight
def find_moves(pos):
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
    queue = deque([('w', 'B3', 'C1', 'A2', 'C3', ())])
    visited = set()
    
    while queue:
        turn, w1, w2, b1, b2, path = queue.popleft()
        
        # Check if the goal is reached
        if (w1, w2) == ('A2', 'C3') and (b1, b2) == ('B3', 'C1'):
            return list(path)
        
        # Generate next moves
        if turn == 'w':
            for move in find_moves(w1):
                new_path = path + (f"w,{w1},{move}",)
                new_state = ('B', move, w2, b1, b2, new_path)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append(new_state)
            for move in find_moves(w2):
                new_path = path + (f"w,{w2},{move}",)
                new_state = ('B', w1, move, b1, b2, new_path)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append(new_state)
        else:
            for move in find_moves(b1):
                new_path = path + (f"B,{b1},{move}",)
                new_state = ('w', w1, w2, move, b2, new_path)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append(new_state)
            for move in find_moves(b2):
                new_path = path + (f"B,{b2},{move}",)
                new_state = ('w', w1, w2, b1, move, new_path)
                if new_state not in visited:
                    visited.add(new_state)
                    queue.append(new_state)
    
    return "No"

# Execute the BFS to find the solution
solution = bfs()
print(solution)