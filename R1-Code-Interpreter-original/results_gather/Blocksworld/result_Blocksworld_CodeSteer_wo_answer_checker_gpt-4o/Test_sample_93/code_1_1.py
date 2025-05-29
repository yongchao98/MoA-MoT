def move_block(state, from_stack, to_stack):
    # Move the top block from one stack to another
    block = state[from_stack].pop()
    state[to_stack].append(block)
    return block

def is_goal_state(state, goal_state):
    # Check if the current state matches the goal state
    return state == goal_state

def dfs(state, goal_state, moves):
    if is_goal_state(state, goal_state):
        return True

    for from_stack in range(len(state)):
        if not state[from_stack]:  # Skip empty stacks
            continue
        for to_stack in range(len(state)):
            if from_stack == to_stack:
                continue
            # Make a move
            block = move_block(state, from_stack, to_stack)
            moves.append(f"Move {block} from {from_stack + 1} to {to_stack + 1}")
            # Recurse
            if dfs(state, goal_state, moves):
                return True
            # Backtrack
            move_block(state, to_stack, from_stack)
            moves.pop()

    return False

# Initial and goal states
initial_state = [['D', 'I'], ['H', 'A', 'C'], ['F', 'B', 'E'], ['G']]
goal_state = [['A', 'B', 'C', 'D', 'F'], ['E', 'H', 'I'], ['G'], []]

# Perform DFS to find the solution
moves = []
dfs(initial_state, goal_state, moves)

# Print the moves
print("<<<" + "\n".join(moves) + ">>>")