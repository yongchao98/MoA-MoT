from collections import deque

# Define the board and initial positions
initial_board = {
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
def find_moves(pos, board):
    x, y = pos_to_coord(pos)
    moves = []
    for dx, dy in knight_moves:
        nx, ny = x + dx, y + dy
        if is_within_board(nx, ny):
            new_pos = coord_to_pos(nx, ny)
            if board[new_pos] is None:
                moves.append(new_pos)
    return moves

# Heuristic: Calculate the distance to the target positions
def heuristic(w1, w2, b1, b2):
    target_w1, target_w2 = 'A2', 'C3'
    target_b1, target_b2 = 'B3', 'C1'
    return (abs(pos_to_coord(w1)[0] - pos_to_coord(target_w1)[0]) +
            abs(pos_to_coord(w1)[1] - pos_to_coord(target_w1)[1]) +
            abs(pos_to_coord(w2)[0] - pos_to_coord(target_w2)[0]) +
            abs(pos_to_coord(w2)[1] - pos_to_coord(target_w2)[1]) +
            abs(pos_to_coord(b1)[0] - pos_to_coord(target_b1)[0]) +
            abs(pos_to_coord(b1)[1] - pos_to_coord(target_b1)[1]) +
            abs(pos_to_coord(b2)[0] - pos_to_coord(target_b2)[0]) +
            abs(pos_to_coord(b2)[1] - pos_to_coord(target_b2)[1]))

# A* search to find the sequence of moves
def a_star_search():
    start = ('w', 'B3', 'C1', 'A2', 'C3', ())
    queue = deque([start])
    visited = set()
    visited.add(start[1:5])
    
    while queue:
        queue = deque(sorted(queue, key=lambda x: heuristic(*x[1:5])))
        turn, w1, w2, b1, b2, path = queue.popleft()
        
        # Check if the goal is reached
        if (w1, w2) == ('A2', 'C3') and (b1, b2) == ('B3', 'C1'):
            return list(path)
        
        # Generate next moves
        current_board = initial_board.copy()
        current_board[w1], current_board[w2] = 'w', 'w'
        current_board[b1], current_board[b2] = 'B', 'B'
        
        if turn == 'w':
            for move in find_moves(w1, current_board):
                new_path = path + (f"w,{w1},{move}",)
                new_state = ('B', move, w2, b1, b2, new_path)
                if new_state[1:5] not in visited:
                    visited.add(new_state[1:5])
                    queue.append(new_state)
            for move in find_moves(w2, current_board):
                new_path = path + (f"w,{w2},{move}",)
                new_state = ('B', w1, move, b1, b2, new_path)
                if new_state[1:5] not in visited:
                    visited.add(new_state[1:5])
                    queue.append(new_state)
        else:
            for move in find_moves(b1, current_board):
                new_path = path + (f"B,{b1},{move}",)
                new_state = ('w', w1, w2, move, b2, new_path)
                if new_state[1:5] not in visited:
                    visited.add(new_state[1:5])
                    queue.append(new_state)
            for move in find_moves(b2, current_board):
                new_path = path + (f"B,{b2},{move}",)
                new_state = ('w', w1, w2, b1, move, new_path)
                if new_state[1:5] not in visited:
                    visited.add(new_state[1:5])
                    queue.append(new_state)
    
    return "No"

# Execute the A* search to find the solution
solution = a_star_search()
print(solution)