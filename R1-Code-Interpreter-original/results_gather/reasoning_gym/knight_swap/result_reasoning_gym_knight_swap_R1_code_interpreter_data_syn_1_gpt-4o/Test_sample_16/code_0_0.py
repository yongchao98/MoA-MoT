from collections import deque

# Define the initial and target positions
initial_positions = {'w': [(0, 2), (2, 0)], 'B': [(3, 1), (3, 0)]}
target_positions = {'w': [(3, 1), (3, 0)], 'B': [(0, 2), (2, 0)]}

# Define the board size
board_size = (4, 3)

# Define knight moves (L-shape)
knight_moves = [(2, 1), (2, -1), (-2, 1), (-2, -1), (1, 2), (1, -2), (-1, 2), (-1, -2)]

# Function to check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < board_size[0] and 0 <= y < board_size[1]

# Function to generate possible moves for a knight
def generate_moves(position):
    x, y = position
    moves = []
    for dx, dy in knight_moves:
        new_x, new_y = x + dx, y + dy
        if is_within_board(new_x, new_y):
            moves.append((new_x, new_y))
    return moves

# BFS to find the sequence of moves
def bfs_solve():
    queue = deque([(initial_positions, [])])
    visited = set()
    visited.add(tuple(sorted(initial_positions['w']) + sorted(initial_positions['B'])))
    
    while queue:
        current_positions, moves = queue.popleft()
        
        # Check if we reached the target positions
        if (sorted(current_positions['w']) == sorted(target_positions['w']) and
            sorted(current_positions['B']) == sorted(target_positions['B'])):
            return moves
        
        # Generate next moves for white knights
        for i, pos in enumerate(current_positions['w']):
            for new_pos in generate_moves(pos):
                if new_pos not in current_positions['w'] and new_pos not in current_positions['B']:
                    new_positions = {'w': current_positions['w'][:], 'B': current_positions['B'][:]}
                    new_positions['w'][i] = new_pos
                    new_positions_tuple = tuple(sorted(new_positions['w']) + sorted(new_positions['B']))
                    if new_positions_tuple not in visited:
                        visited.add(new_positions_tuple)
                        queue.append((new_positions, moves + [f"w,{pos},{new_pos}"]))
        
        # Generate next moves for black knights
        for i, pos in enumerate(current_positions['B']):
            for new_pos in generate_moves(pos):
                if new_pos not in current_positions['w'] and new_pos not in current_positions['B']:
                    new_positions = {'w': current_positions['w'][:], 'B': current_positions['B'][:]}
                    new_positions['B'][i] = new_pos
                    new_positions_tuple = tuple(sorted(new_positions['w']) + sorted(new_positions['B']))
                    if new_positions_tuple not in visited:
                        visited.add(new_positions_tuple)
                        queue.append((new_positions, moves + [f"B,{pos},{new_pos}"]))
    
    return "No"

# Find the solution
solution = bfs_solve()
print(solution)