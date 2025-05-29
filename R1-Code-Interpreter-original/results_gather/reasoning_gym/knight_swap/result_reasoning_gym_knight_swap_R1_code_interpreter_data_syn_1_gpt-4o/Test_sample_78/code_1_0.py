from collections import deque

# Define the board and initial positions
board = {
    'A3': 'B', 'B1': 'B',
    'C2': 'w', 'C1': 'w'
}

# Define the target positions
target_positions = {
    'w': ['A3', 'B1'],
    'B': ['C2', 'C1']
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

# Check if a move is valid
def is_valid_move(pos, move, current_board):
    x, y = pos_to_coord(pos)
    dx, dy = move
    new_x, new_y = x + dx, y + dy
    new_pos = coord_to_pos((new_x, new_y))
    return 0 <= new_x < 4 and 0 <= new_y < 3 and new_pos not in current_board

# BFS to find the sequence of moves
def bfs():
    queue = deque([('w', board, [])])
    visited = set()
    
    while queue:
        color, current_board, moves = queue.popleft()
        
        # Check if the current board matches the target positions
        if all(pos in current_board and current_board[pos] == color for pos in target_positions[color]):
            return moves
        
        # Alternate color
        next_color = 'B' if color == 'w' else 'w'
        
        # Try all possible moves for the current color
        for pos, knight in current_board.items():
            if knight == color:
                for move in knight_moves:
                    if is_valid_move(pos, move, current_board):
                        new_pos = coord_to_pos((pos_to_coord(pos)[0] + move[0], pos_to_coord(pos)[1] + move[1]))
                        new_board = current_board.copy()
                        new_board[new_pos] = color
                        del new_board[pos]
                        new_moves = moves + [f"{color},{pos},{new_pos}"]
                        board_tuple = tuple(sorted(new_board.items()))
                        if board_tuple not in visited:
                            visited.add(board_tuple)
                            queue.append((next_color, new_board, new_moves))
    
    return "No"

# Find the solution
solution = bfs()
print(solution)