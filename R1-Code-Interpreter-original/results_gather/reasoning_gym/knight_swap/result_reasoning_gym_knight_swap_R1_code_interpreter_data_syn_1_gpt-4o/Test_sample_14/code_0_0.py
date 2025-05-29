from collections import deque

# Define the board and initial positions
board = {
    'A3': '.', 'B3': 'B', 'C3': 'B', 'D3': 'w',
    'A2': '.', 'B2': 'B', 'C2': '.', 'D2': '.',
    'A1': '.', 'B1': 'w', 'C1': '.', 'D1': '.'
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
def coord_to_pos(coord):
    col, row = coord
    return chr(col + ord('A')) + str(row + 1)

# Check if a move is within the board
def is_within_board(coord):
    col, row = coord
    return 0 <= col < 4 and 0 <= row < 3

# Get possible moves for a knight
def get_possible_moves(pos):
    col, row = pos_to_coord(pos)
    possible_moves = []
    for move in knight_moves:
        new_pos = (col + move[0], row + move[1])
        if is_within_board(new_pos):
            new_pos_str = coord_to_pos(new_pos)
            if board[new_pos_str] == '.':
                possible_moves.append(new_pos_str)
    return possible_moves

# BFS to find the sequence of moves
def bfs():
    queue = deque([('w', 'D3', 'B1', [])])  # Start with white knight at D3
    visited = set()
    
    while queue:
        color, pos1, pos2, moves = queue.popleft()
        
        # Check if the current state is the target state
        if (pos1, pos2) == ('C3', 'B2'):
            return moves
        
        # Alternate color
        next_color = 'B' if color == 'w' else 'w'
        
        # Get possible moves for the current knight
        for move in get_possible_moves(pos1):
            new_moves = moves + [f"{color},{pos1},{move}"]
            new_state = (next_color, move, pos2) if color == 'w' else (next_color, pos1, move)
            if new_state not in visited:
                visited.add(new_state)
                queue.append((*new_state, new_moves))
    
    return "No"

# Run the BFS to find the solution
result = bfs()
print(result)