from collections import deque

# Define the board and initial positions
initial_positions = {
    'w': [(2, 0), (3, 2)],  # C1, D3
    'B': [(1, 2), (2, 1)]   # B3, C2
}

# Target positions for the swap
target_positions = {
    'w': [(1, 2), (2, 1)],  # B3, C2
    'B': [(2, 0), (3, 2)]   # C1, D3
}

# Possible knight moves (L-shape)
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert board positions to coordinates
def pos_to_coord(pos):
    col, row = pos
    return (ord(col) - ord('A'), int(row) - 1)

# Convert coordinates to board positions
def coord_to_pos(coord):
    x, y = coord
    return f"{chr(x + ord('A'))}{y + 1}"

# Check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 4 and 0 <= y < 3

# BFS to find the sequence of moves
def bfs_knight_swap():
    queue = deque()
    visited = set()
    
    # Initial state
    initial_state = (tuple(initial_positions['w']), tuple(initial_positions['B']))
    queue.append((initial_state, []))
    visited.add(initial_state)
    
    while queue:
        (w_pos, b_pos), moves = queue.popleft()
        
        # Check if we reached the target state
        if set(w_pos) == set(target_positions['w']) and set(b_pos) == set(target_positions['B']):
            return moves
        
        # Generate next states for white knights
        for i, (x, y) in enumerate(w_pos):
            for dx, dy in knight_moves:
                nx, ny = x + dx, y + dy
                if is_within_board(nx, ny) and (nx, ny) not in w_pos and (nx, ny) not in b_pos:
                    new_w_pos = list(w_pos)
                    new_w_pos[i] = (nx, ny)
                    new_state = (tuple(new_w_pos), b_pos)
                    if new_state not in visited:
                        visited.add(new_state)
                        move = f"w,{coord_to_pos((x, y))},{coord_to_pos((nx, ny))}"
                        queue.append((new_state, moves + [move]))
        
        # Generate next states for black knights
        for i, (x, y) in enumerate(b_pos):
            for dx, dy in knight_moves:
                nx, ny = x + dx, y + dy
                if is_within_board(nx, ny) and (nx, ny) not in w_pos and (nx, ny) not in b_pos:
                    new_b_pos = list(b_pos)
                    new_b_pos[i] = (nx, ny)
                    new_state = (w_pos, tuple(new_b_pos))
                    if new_state not in visited:
                        visited.add(new_state)
                        move = f"B,{coord_to_pos((x, y))},{coord_to_pos((nx, ny))}"
                        queue.append((new_state, moves + [move]))
    
    return "No"

# Find the sequence of moves
result = bfs_knight_swap()
print(result)