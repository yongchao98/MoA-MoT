from collections import deque

# Define the board and initial positions
initial_positions = {
    'w': [(2, 0), (0, 2)],  # C3, A2
    'B': [(0, 0), (2, 2)]   # A1, C1
}

# Target positions
target_positions = {
    'w': [(0, 0), (2, 2)],  # A1, C1
    'B': [(2, 0), (0, 2)]   # C3, A2
}

# Possible knight moves (L-shape)
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert position to board notation
def pos_to_notation(pos):
    col = chr(ord('A') + pos[1])
    row = str(pos[0] + 1)
    return col + row

# Check if a position is within the board
def is_within_board(pos):
    return 0 <= pos[0] < 3 and 0 <= pos[1] < 4

# BFS to find the sequence of moves
def find_moves():
    queue = deque()
    visited = set()
    
    # Initial state
    initial_state = (tuple(initial_positions['w']), tuple(initial_positions['B']), 'w')
    queue.append((initial_state, []))
    visited.add(initial_state)
    
    while queue:
        (w_pos, b_pos, turn), moves = queue.popleft()
        
        # Check if we reached the target
        if (w_pos == tuple(target_positions['w']) and
            b_pos == tuple(target_positions['B'])):
            return moves
        
        # Determine current positions and next turn
        current_positions = w_pos if turn == 'w' else b_pos
        next_turn = 'B' if turn == 'w' else 'w'
        
        # Try all possible moves for the current knight
        for i, knight_pos in enumerate(current_positions):
            for move in knight_moves:
                new_pos = (knight_pos[0] + move[0], knight_pos[1] + move[1])
                
                if is_within_board(new_pos) and new_pos not in w_pos and new_pos not in b_pos:
                    # Create new positions
                    new_positions = list(current_positions)
                    new_positions[i] = new_pos
                    new_positions = tuple(new_positions)
                    
                    # Create new state
                    if turn == 'w':
                        new_state = (new_positions, b_pos, next_turn)
                    else:
                        new_state = (w_pos, new_positions, next_turn)
                    
                    if new_state not in visited:
                        visited.add(new_state)
                        new_move = f"{turn},{pos_to_notation(knight_pos)},{pos_to_notation(new_pos)}"
                        queue.append((new_state, moves + [new_move]))
    
    return "No"

# Find the sequence of moves
result = find_moves()
print(result)