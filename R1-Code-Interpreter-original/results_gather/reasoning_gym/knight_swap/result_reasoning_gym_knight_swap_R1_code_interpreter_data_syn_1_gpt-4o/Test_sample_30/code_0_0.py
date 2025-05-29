from collections import deque

# Define the initial and target positions
initial_positions = {'w': [(0, 1), (1, 1)], 'B': [(2, 0), (3, 0)]}
target_positions = {'w': [(2, 0), (3, 0)], 'B': [(0, 1), (1, 1)]}

# Define the board size
board_size = (4, 3)

# Define knight moves (L-shape)
knight_moves = [(2, 1), (2, -1), (-2, 1), (-2, -1), (1, 2), (1, -2), (-1, 2), (-1, -2)]

# Function to check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < board_size[0] and 0 <= y < board_size[1]

# Function to get the board state as a tuple
def get_board_state(positions):
    board = [['.' for _ in range(board_size[0])] for _ in range(board_size[1])]
    for color, pos_list in positions.items():
        for x, y in pos_list:
            board[y][x] = color
    return tuple(tuple(row) for row in board)

# BFS to find the sequence of moves
def bfs_knight_swap():
    # Initial state
    initial_state = get_board_state(initial_positions)
    target_state = get_board_state(target_positions)
    
    # Queue for BFS
    queue = deque([(initial_positions, [])])
    visited = set()
    visited.add(initial_state)
    
    while queue:
        current_positions, moves = queue.popleft()
        
        # Check if we reached the target state
        if get_board_state(current_positions) == target_state:
            return moves
        
        # Get all possible moves for the current player
        for color in ['w', 'B']:
            for i, (x, y) in enumerate(current_positions[color]):
                for dx, dy in knight_moves:
                    nx, ny = x + dx, y + dy
                    if is_within_board(nx, ny):
                        # Check if the target square is empty
                        if all((nx, ny) not in pos_list for pos_list in current_positions.values()):
                            # Make the move
                            new_positions = {c: list(pos_list) for c, pos_list in current_positions.items()}
                            new_positions[color][i] = (nx, ny)
                            new_state = get_board_state(new_positions)
                            
                            if new_state not in visited:
                                visited.add(new_state)
                                new_moves = moves + [f"{color},{chr(65+x)}{3-y},{chr(65+nx)}{3-ny}"]
                                queue.append((new_positions, new_moves))
    
    return "No"

# Run the BFS to find the solution
result = bfs_knight_swap()
print(result)