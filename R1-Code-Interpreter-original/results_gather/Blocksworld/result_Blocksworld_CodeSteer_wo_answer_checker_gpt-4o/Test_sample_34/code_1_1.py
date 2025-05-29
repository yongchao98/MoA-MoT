def move_block(state, from_stack, to_stack):
    """Move the top block from one stack to another."""
    block = state[from_stack].pop()
    state[to_stack].append(block)
    return block

def is_goal_state(state, goal_state):
    """Check if the current state matches the goal state."""
    return state == goal_state

def solve_blocks(state, goal_state, moves, visited):
    """Recursive function to solve the blocksworld problem."""
    if is_goal_state(state, goal_state):
        return True

    state_tuple = tuple(tuple(stack) for stack in state)
    if state_tuple in visited:
        return False

    visited.add(state_tuple)

    for from_stack in range(len(state)):
        if not state[from_stack]:
            continue
        for to_stack in range(len(state)):
            if from_stack == to_stack:
                continue
            # Make a move
            block = move_block(state, from_stack, to_stack)
            moves.append(f"Move {block} from {from_stack + 1} to {to_stack + 1}")

            # Recurse
            if solve_blocks(state, goal_state, moves, visited):
                return True

            # Backtrack
            moves.pop()
            move_block(state, to_stack, from_stack)

    return False

# Initial and goal states
initial_state = [['G', 'E', 'B'], ['A'], ['C', 'D', 'H', 'F']]
goal_state = [['D', 'F', 'G'], ['A', 'C', 'E', 'H'], ['B']]

# Solve the problem
moves = []
visited = set()
if solve_blocks(initial_state, goal_state, moves, visited):
    print("<<<")
    for move in moves:
        print(move)
    print(">>>")
else:
    print("No solution found.")