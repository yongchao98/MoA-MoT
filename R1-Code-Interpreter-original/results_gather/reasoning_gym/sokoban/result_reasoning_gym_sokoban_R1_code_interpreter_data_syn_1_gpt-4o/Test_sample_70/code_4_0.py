def is_deadlock(state):
    # Check for simple deadlocks where a box is in a corner not on a goal
    for box in state['boxes']:
        if box not in state['goals']:
            x, y = box
            if ((x - 1, y) in state['walls'] or (x + 1, y) in state['walls']) and \
               ((x, y - 1) in state['walls'] or (x, y + 1) in state['walls']):
                return True
    return False

def dfs(state, depth, path, visited):
    if is_goal_state(state):
        return path
    if depth == 0 or is_deadlock(state):
        return None
    
    for direction in moves:
        new_state = move_player(state, direction)
        if new_state is None:
            continue
        
        state_signature = (new_state['player'], frozenset(new_state['boxes']))
        if state_signature not in visited:
            visited.add(state_signature)
            result = dfs(new_state, depth - 1, path + direction, visited)
            if result is not None:
                return result
            visited.remove(state_signature)
    
    return None

def solve_sokoban_iterative_deepening(initial_state):
    depth = 0
    while True:
        visited = set()
        visited.add((initial_state['player'], frozenset(initial_state['boxes'])))
        result = dfs(initial_state, depth, "", visited)
        if result is not None:
            return result
        depth += 1

solution = solve_sokoban_iterative_deepening(initial_state)
print(solution)