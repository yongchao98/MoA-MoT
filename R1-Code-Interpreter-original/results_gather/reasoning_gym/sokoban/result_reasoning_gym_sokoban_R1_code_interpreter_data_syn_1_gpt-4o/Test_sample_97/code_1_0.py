def depth_limited_search(state, depth, path, visited):
    if is_solved(state[1]):
        return path
    if depth == 0:
        return None
    visited.add(state)
    for move_name, move in moves.items():
        if is_valid_move(state[0], state[1], move):
            new_state = apply_move(state[0], state[1], move)
            if new_state not in visited and not is_deadlock(new_state[1]):
                result = depth_limited_search(new_state, depth - 1, path + move_name, visited)
                if result is not None:
                    return result
    visited.remove(state)
    return None

def iterative_deepening_search(initial_state):
    depth = 0
    while True:
        visited = set()
        result = depth_limited_search(initial_state, depth, "", visited)
        if result is not None:
            return result
        depth += 1

# Define the initial state of the puzzle
initial_state = (
    (8, 2),  # Player position
    frozenset([(2, 2), (3, 3), (3, 4), (5, 2), (7, 2), (7, 3), (7, 7)])  # Box positions
)

solution = iterative_deepening_search(initial_state)
print(solution)