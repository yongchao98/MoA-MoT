def move_block(state, from_stack, to_stack):
    """Move the top block from one stack to another."""
    block = state[from_stack].pop()
    state[to_stack].append(block)
    return block

def is_goal_state(state, goal_state):
    """Check if the current state matches the goal state."""
    return state == goal_state

def dfs(state, goal_state, path):
    """Perform a depth-first search to find a sequence of moves."""
    if is_goal_state(state, goal_state):
        return path

    for from_stack in range(len(state)):
        if not state[from_stack]:
            continue
        for to_stack in range(len(state)):
            if from_stack != to_stack:
                # Copy the current state to explore this move
                new_state = [stack[:] for stack in state]
                block = move_block(new_state, from_stack, to_stack)
                new_path = path + [f"Move {block} from stack{from_stack + 1} to stack{to_stack + 1}"]
                result = dfs(new_state, goal_state, new_path)
                if result:
                    return result
    return None

# Initial and goal states
initial_state = [['B', 'G', 'I'], ['F', 'H', 'D'], ['E', 'C'], ['A']]
goal_state = [['B', 'D', 'E', 'I'], ['F', 'G', 'H'], ['A', 'C']]

# Find the sequence of moves
solution = dfs(initial_state, goal_state, [])
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")