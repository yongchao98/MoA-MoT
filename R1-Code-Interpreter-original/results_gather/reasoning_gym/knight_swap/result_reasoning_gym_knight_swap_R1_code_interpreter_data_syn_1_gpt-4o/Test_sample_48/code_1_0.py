from collections import deque

# Define the initial and target positions
initial_positions = {'w': [(2, 2), (1, 1)], 'B': [(1, 2), (0, 1)]}
target_positions = {'w': [(1, 2), (0, 1)], 'B': [(2, 2), (1, 1)]}

# Define the board size
board_size = (4, 3)

# Define knight moves
knight_moves = [(2, 1), (2, -1), (-2, 1), (-2, -1), (1, 2), (1, -2), (-1, 2), (-1, -2)]

# Function to check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < board_size[0] and 0 <= y < board_size[1]

# Function to get the next positions for a knight
def get_knight_moves(pos, occupied):
    x, y = pos
    moves = []
    for dx, dy in knight_moves:
        nx, ny = x + dx, y + dy
        if is_within_board(nx, ny) and (nx, ny) not in occupied:
            moves.append((nx, ny))
    return moves

# BFS to find the sequence of moves
def bfs_knight_swap():
    queue = deque([(initial_positions, [])])
    visited = set()
    visited.add((tuple(initial_positions['w']), tuple(initial_positions['B'])))
    
    while queue:
        current_positions, moves = queue.popleft()
        
        # Check if we reached the target positions
        if current_positions == target_positions:
            return moves
        
        # Get the current positions of white and black knights
        white_positions = current_positions['w']
        black_positions = current_positions['B']
        
        # Get all occupied positions
        occupied = set(white_positions + black_positions)
        
        # Try all possible moves for white knights
        for i, pos in enumerate(white_positions):
            for new_pos in get_knight_moves(pos, occupied):
                new_white_positions = list(white_positions)
                new_white_positions[i] = new_pos
                new_occupied = set(new_white_positions + list(black_positions))
                
                # Try all possible moves for black knights
                for j, b_pos in enumerate(black_positions):
                    for new_b_pos in get_knight_moves(b_pos, new_occupied):
                        new_black_positions = list(black_positions)
                        new_black_positions[j] = new_b_pos
                        
                        new_state = {'w': tuple(new_white_positions), 'B': tuple(new_black_positions)}
                        if (new_state['w'], new_state['B']) not in visited:
                            visited.add((new_state['w'], new_state['B']))
                            new_moves = moves + [f"w,{pos},{new_pos}", f"B,{b_pos},{new_b_pos}"]
                            queue.append((new_state, new_moves))
    
    return "No"

# Execute the BFS to find the solution
result = bfs_knight_swap()
print(result)