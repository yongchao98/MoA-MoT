from collections import deque

# Define the board positions
initial_positions = {
    'w': [(1, 1), (2, 1)],  # B1, C1
    'B': [(2, 3), (0, 0)]   # C3, A1
}

# Target positions
target_positions = {
    'w': [(2, 3), (0, 0)],  # C3, A1
    'B': [(1, 1), (2, 1)]   # B1, C1
}

# Possible knight moves (L-shape)
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert board coordinates to chess notation
def to_chess_notation(x, y):
    return chr(ord('A') + x) + str(y + 1)

# Check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 4 and 0 <= y < 4

# BFS to find the sequence of moves
def bfs_knight_swap():
    queue = deque()
    visited = set()
    
    # Initial state
    initial_state = (tuple(initial_positions['w']), tuple(initial_positions['B']))
    queue.append((initial_state, []))
    visited.add(initial_state)
    
    while queue:
        (w_positions, b_positions), moves = queue.popleft()
        
        # Check if we reached the target state
        if set(w_positions) == set(target_positions['w']) and set(b_positions) == set(target_positions['B']):
            return moves
        
        # Generate next states for white knights
        for i, (wx, wy) in enumerate(w_positions):
            for dx, dy in knight_moves:
                nx, ny = wx + dx, wy + dy
                if is_within_board(nx, ny) and (nx, ny) not in w_positions and (nx, ny) not in b_positions:
                    new_w_positions = list(w_positions)
                    new_w_positions[i] = (nx, ny)
                    new_state = (tuple(new_w_positions), b_positions)
                    if new_state not in visited:
                        visited.add(new_state)
                        move = f"w,{to_chess_notation(wx, wy)},{to_chess_notation(nx, ny)}"
                        queue.append((new_state, moves + [move]))
        
        # Generate next states for black knights
        for i, (bx, by) in enumerate(b_positions):
            for dx, dy in knight_moves:
                nx, ny = bx + dx, by + dy
                if is_within_board(nx, ny) and (nx, ny) not in w_positions and (nx, ny) not in b_positions:
                    new_b_positions = list(b_positions)
                    new_b_positions[i] = (nx, ny)
                    new_state = (w_positions, tuple(new_b_positions))
                    if new_state not in visited:
                        visited.add(new_state)
                        move = f"B,{to_chess_notation(bx, by)},{to_chess_notation(nx, ny)}"
                        queue.append((new_state, moves + [move]))
    
    return "No"

# Find the solution
solution = bfs_knight_swap()
print(solution)