def move_block(state, from_stack, to_stack):
    """Move the top block from one stack to another."""
    new_state = [list(stack) for stack in state]  # Deep copy of the state
    block = new_state[from_stack].pop()
    new_state[to_stack].append(block)
    return new_state

def is_goal(state, goal):
    """Check if the current state matches the goal state."""
    return state == goal

def dfs(state, goal, path, visited):
    """Depth-First Search to find the sequence of moves to reach the goal."""
    if is_goal(state, goal):
        return path

    visited.add(tuple(tuple(stack) for stack in state))  # Mark the state as visited

    for from_stack in range(3):
        if not state[from_stack]:  # Skip empty stacks
            continue
        for to_stack in range(3):
            if from_stack == to_stack:
                continue
            new_state = move_block(state, from_stack, to_stack)
            new_state_tuple = tuple(tuple(stack) for stack in new_state)
            if new_state_tuple not in visited:
                result = dfs(new_state, goal, path + [f"Move {state[from_stack][-1]} from stack{from_stack+1} to stack{to_stack+1}"], visited)
                if result:
                    return result
    return None

# Initial and goal states
initial_state = [['C', 'D', 'B'], ['E'], ['F', 'A']]
goal_state = [['A', 'C', 'F'], ['B'], ['D', 'E']]

# Perform DFS to find the solution
solution = dfs(initial_state, goal_state, [], set())

# Print the solution
if solution:
    print("<<<")
    for move in solution:
        print(move)
    print(">>>")
else:
    print("No solution found.")