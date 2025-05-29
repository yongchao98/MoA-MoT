def move_block(state, from_stack, to_stack):
    """Move the top block from one stack to another."""
    block = state[from_stack].pop()
    state[to_stack].append(block)
    return block

def is_goal_state(state, goal_state):
    """Check if the current state matches the goal state."""
    return state == goal_state

def state_to_tuple(state):
    """Convert the state to a tuple for hashing."""
    return tuple(tuple(stack) for stack in state)

def iddfs(state, goal_state, max_depth):
    """Perform an iterative deepening depth-first search."""
    for depth in range(max_depth):
        visited = set()
        path = []
        result = dfs(state, goal_state, path, visited, depth)
        if result:
            return result
    return None

def dfs(state, goal_state, path, visited, depth):
    """Perform a depth-first search with a depth limit."""
    if depth == 0:
        return None
    if is_goal_state(state, goal_state):
        return path

    state_tuple = state_to_tuple(state)
    if state_tuple in visited:
        return None
    visited.add(state_tuple)

    for from_stack in range(len(state)):
        if not state[from_stack]:
            continue
        for to_stack in range(len(state)):
            if from_stack != to_stack:
                # Copy the current state to explore this move
                new_state = [stack[:] for stack in state]
                block = move_block(new_state, from_stack, to_stack)
                new_path = path + [f"Move {block} from stack{from_stack + 1} to stack{to_stack + 1}"]
                result = dfs(new_state, goal_state, new_path, visited, depth - 1)
                if result:
                    return result
    return None

# Initial and goal states
initial_state = [['B', 'G', 'I'], ['F', 'H', 'D'], ['E', 'C'], ['A']]
goal_state = [['B', 'D', 'E', 'I'], ['F', 'G', 'H'], ['A', 'C']]

# Find the sequence of moves with a reasonable depth limit
solution = iddfs(initial_state, goal_state, max_depth=20)
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")