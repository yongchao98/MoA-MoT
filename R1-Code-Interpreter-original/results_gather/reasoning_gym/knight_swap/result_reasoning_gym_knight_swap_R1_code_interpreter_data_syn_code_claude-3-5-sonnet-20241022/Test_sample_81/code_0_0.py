from collections import deque

def get_knight_moves(pos, size=4):
    moves = [
        (2,1), (2,-1), (-2,1), (-2,-1),
        (1,2), (1,-2), (-1,2), (-1,-2)
    ]
    x, y = pos
    valid = []
    for dx, dy in moves:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < size and 0 <= new_y < size:
            valid.append((new_x, new_y))
    return valid

def pos_to_coord(pos):
    return f"{chr(65+pos[1])}{pos[0]+1}"

def state_to_grid(state):
    return [list(state[i:i+4]) for i in range(0, 16, 4)]

def get_pieces_positions(grid):
    positions = {'w': [], 'B': []}
    for i in range(4):
        for j in range(4):
            if grid[i][j] in 'wB':
                positions[grid[i][j]].append((i, j))
    return positions

def is_target_reached(positions):
    # Check if white knights are at (2,2), (1,3) and black at (2,3), (1,2)
    w_target = {(2,2), (1,3)}
    b_target = {(2,3), (1,2)}
    return set(positions['w']) == b_target and set(positions['B']) == w_target

def make_move(grid, from_pos, to_pos, piece):
    new_grid = [row[:] for row in grid]
    new_grid[from_pos[0]][from_pos[1]] = '.'
    new_grid[to_pos[0]][to_pos[1]] = piece
    return [''.join(row) for row in new_grid]

def solve_knight_swap():
    initial_state = "...B" + "..wB" + "..." + "...."
    queue = deque([(initial_state, [], True)])  # state, moves, is_white_turn
    visited = {initial_state}
    
    while queue:
        state, moves, is_white_turn = queue.popleft()
        grid = state_to_grid(state)
        positions = get_pieces_positions(grid)
        
        if is_target_reached(positions):
            return moves
        
        current_piece = 'w' if is_white_turn else 'B'
        for piece_pos in positions[current_piece]:
            for new_pos in get_knight_moves(piece_pos):
                if grid[new_pos[0]][new_pos[1]] == '.':
                    new_grid = make_move(grid, piece_pos, new_pos, current_piece)
                    new_state = ''.join([''.join(row) for row in new_grid])
                    
                    if new_state not in visited:
                        visited.add(new_state)
                        new_move = f"{current_piece},{pos_to_coord(piece_pos)},{pos_to_coord(new_pos)}"
                        queue.append((new_state, moves + [new_move], not is_white_turn))
    
    return None

solution = solve_knight_swap()
if solution is None:
    print("<<<No>>>")
else:
    print(f"<<<{solution}>>>")