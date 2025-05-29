from heapq import heappush, heappop
import copy

def get_initial_state():
    puzzle = [
        "+ + + + + + + + + +",
        "+ - - X - - - @ * +",
        "+ - - - @ X @ @ @ +",
        "+ - - - X + X $ - +",
        "+ + - - + $ + + X +",
        "+ X - - + + + + $ +",
        "+ - @ - + + + + + +",
        "+ - - - - + + + + +",
        "+ - - - - + + + $ +",
        "+ + + + + + + + + +"
    ]
    return [list(row.replace(" ", "")) for row in puzzle]

def manhattan_distance(pos1, pos2):
    return abs(pos1[0] - pos2[0]) + abs(pos1[1] - pos2[1])

def get_positions(state, chars):
    positions = []
    for i in range(len(state)):
        for j in range(len(state[i])):
            if state[i][j] in chars:
                positions.append((i, j))
    return positions

def calculate_heuristic(state):
    # Get positions of boxes and goals
    boxes = get_positions(state, ['@', '$'])
    goals = get_positions(state, ['X', '$', '*'])
    
    # Calculate minimum sum of distances from each box to a goal
    total_distance = 0
    used_goals = set()
    
    for box in boxes:
        min_dist = float('inf')
        best_goal = None
        for goal in goals:
            if goal not in used_goals:
                dist = manhattan_distance(box, goal)
                if dist < min_dist:
                    min_dist = dist
                    best_goal = goal
        if best_goal:
            used_goals.add(best_goal)
            total_distance += min_dist
    
    return total_distance

def get_player_pos(state):
    for i in range(len(state)):
        for j in range(len(state[i])):
            if state[i][j] in ['@', '*']:
                return (i, j)
    return None

def get_valid_moves(state):
    moves = []
    player_pos = get_player_pos(state)
    if not player_pos:
        return moves
    
    directions = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}
    
    for direction, (dr, dc) in directions.items():
        new_state = [row[:] for row in state]
        r, c = player_pos
        new_r, new_c = r + dr, c + dc
        
        if 0 <= new_r < len(state) and 0 <= new_c < len(state[0]):
            curr_char = state[r][c]
            next_char = state[new_r][new_c]
            
            # Moving to empty space or goal
            if next_char in ['-', 'X']:
                new_state[r][c] = 'X' if curr_char == '*' else '-'
                new_state[new_r][new_c] = '*' if next_char == 'X' else '@'
                moves.append((direction, new_state))
            
            # Pushing a box
            elif next_char in ['@', '$']:
                push_r, push_c = new_r + dr, new_c + dc
                if 0 <= push_r < len(state) and 0 <= push_c < len(state[0]):
                    push_char = state[push_r][push_c]
                    if push_char in ['-', 'X']:
                        new_state[r][c] = 'X' if curr_char == '*' else '-'
                        new_state[new_r][new_c] = '*' if next_char == '$' else '@'
                        new_state[push_r][push_c] = '$' if push_char == 'X' else '@'
                        moves.append((direction, new_state))
    
    return moves

def is_goal_state(state):
    for row in state:
        for cell in row:
            if cell == '@':  # If there's any box not on a goal
                return False
    return True

def state_to_string(state):
    return ''.join(''.join(row) for row in state)

def solve_sokoban():
    initial_state = get_initial_state()
    initial_h = calculate_heuristic(initial_state)
    
    # Priority queue: (f_score, moves_count, path, state)
    pq = [(initial_h, 0, "", initial_state)]
    visited = {state_to_string(initial_state)}
    max_moves = 50  # Limit search depth
    
    while pq:
        _, moves_count, path, state = heappop(pq)
        
        if moves_count > max_moves:
            continue
            
        if is_goal_state(state):
            return path
        
        for move, new_state in get_valid_moves(state):
            state_str = state_to_string(new_state)
            if state_str not in visited:
                visited.add(state_str)
                h = calculate_heuristic(new_state)
                f = h + moves_count + 1
                heappush(pq, (f, moves_count + 1, path + move, new_state))
    
    return None

solution = solve_sokoban()
print(solution if solution else "No solution found")