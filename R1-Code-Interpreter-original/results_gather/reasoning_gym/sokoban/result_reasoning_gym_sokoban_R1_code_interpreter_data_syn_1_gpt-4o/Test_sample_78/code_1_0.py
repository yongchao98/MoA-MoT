from heapq import heappop, heappush

# Define the initial state of the puzzle
initial_state = (
    (1, 3),  # Player position
    ((2, 3), (2, 6), (2, 7), (4, 6))  # Box positions
)

# Define the goal positions
goal_positions = {(1, 1), (2, 4), (4, 4), (5, 7)}

# Define the walls
walls = {
    (0, 0), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8),
    (1, 0), (1, 8),
    (2, 0), (2, 8),
    (3, 0), (3, 4), (3, 5), (3, 8),
    (4, 0), (4, 1), (4, 2), (4, 3), (4, 8),
    (5, 0), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 8),
    (6, 0), (6, 1), (6, 2), (6, 3), (6, 4), (6, 5), (6, 6), (6, 7), (6, 8)
}

# Define possible moves
moves = {
    'U': (-1, 0),
    'D': (1, 0),
    'L': (0, -1),
    'R': (0, 1)
}

def is_goal_state(box_positions):
    return all(pos in goal_positions for pos in box_positions)

def is_valid_position(position, box_positions):
    return position not in walls and position not in box_positions

def move_player(state, move):
    player_pos, box_positions = state
    dx, dy = moves[move]
    new_player_pos = (player_pos[0] + dx, player_pos[1] + dy)
    
    if not is_valid_position(new_player_pos, box_positions):
        return None
    
    new_box_positions = set(box_positions)
    if new_player_pos in box_positions:
        new_box_pos = (new_player_pos[0] + dx, new_player_pos[1] + dy)
        if not is_valid_position(new_box_pos, box_positions):
            return None
        new_box_positions.remove(new_player_pos)
        new_box_positions.add(new_box_pos)
    
    return (new_player_pos, tuple(sorted(new_box_positions)))

def heuristic(box_positions):
    return sum(min(abs(bx - gx) + abs(by - gy) for gx, gy in goal_positions) for bx, by in box_positions)

def solve_sokoban(initial_state):
    priority_queue = []
    heappush(priority_queue, (0, initial_state, ""))
    visited = set()
    visited.add(initial_state)
    
    while priority_queue:
        _, current_state, path = heappop(priority_queue)
        
        if is_goal_state(current_state[1]):
            return path
        
        for move in moves:
            new_state = move_player(current_state, move)
            if new_state and new_state not in visited:
                visited.add(new_state)
                cost = len(path) + 1 + heuristic(new_state[1])
                heappush(priority_queue, (cost, new_state, path + move))
    
    return None

solution = solve_sokoban(initial_state)
print(solution)