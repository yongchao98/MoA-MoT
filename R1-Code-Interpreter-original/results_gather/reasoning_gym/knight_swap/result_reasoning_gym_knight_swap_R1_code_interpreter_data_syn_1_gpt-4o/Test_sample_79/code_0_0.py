from collections import deque

# Define the initial and target positions
initial_positions = {'w': [(1, 2), (2, 2)], 'B': [(3, 2), (2, 0)]}
target_positions = {'w': [(3, 2), (2, 0)], 'B': [(1, 2), (2, 2)]}

# Define the board size
board_size = (4, 3)

# Define knight moves
knight_moves = [(2, 1), (2, -1), (-2, 1), (-2, -1), (1, 2), (1, -2), (-1, 2), (-1, -2)]

# Function to check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < board_size[0] and 0 <= y < board_size[1]

# Function to generate possible moves for a knight
def generate_moves(position, occupied_positions):
    x, y = position
    possible_moves = []
    for dx, dy in knight_moves:
        new_x, new_y = x + dx, y + dy
        if is_within_board(new_x, new_y) and (new_x, new_y) not in occupied_positions:
            possible_moves.append((new_x, new_y))
    return possible_moves

# BFS to find the sequence of moves
def bfs_knight_swap():
    queue = deque([(initial_positions, [])])
    visited = set()
    visited.add(tuple(sorted(initial_positions['w']) + sorted(initial_positions['B'])))
    
    while queue:
        current_positions, moves = queue.popleft()
        
        # Check if we reached the target positions
        if (sorted(current_positions['w']) == sorted(target_positions['w']) and
            sorted(current_positions['B']) == sorted(target_positions['B'])):
            return moves
        
        # Determine whose turn it is
        turn = 'w' if len(moves) % 2 == 0 else 'B'
        
        # Generate all possible moves for the current turn
        for i, knight_pos in enumerate(current_positions[turn]):
            occupied_positions = set(current_positions['w'] + current_positions['B'])
            possible_moves = generate_moves(knight_pos, occupied_positions)
            
            for new_pos in possible_moves:
                new_positions = {color: list(pos) for color, pos in current_positions.items()}
                new_positions[turn][i] = new_pos
                new_positions_tuple = tuple(sorted(new_positions['w']) + sorted(new_positions['B']))
                
                if new_positions_tuple not in visited:
                    visited.add(new_positions_tuple)
                    new_move = f"{turn},{chr(65 + knight_pos[0])}{3 - knight_pos[1]},{chr(65 + new_pos[0])}{3 - new_pos[1]}"
                    queue.append((new_positions, moves + [new_move]))
    
    return "No"

# Execute the BFS function
result = bfs_knight_swap()
print(result)