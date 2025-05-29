from collections import deque

def move_block(state, from_stack, to_stack):
    # Create a new state with the block moved
    new_state = [list(stack) for stack in state]
    block = new_state[from_stack].pop()
    new_state[to_stack].append(block)
    return tuple(tuple(stack) for stack in new_state), block

def is_goal_state(state, goal_state):
    return state == goal_state

def bfs(initial_state, goal_state):
    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(initial_state)

    while queue:
        current_state, moves = queue.popleft()

        if is_goal_state(current_state, goal_state):
            return moves

        for from_stack in range(len(current_state)):
            if not current_state[from_stack]:  # Skip empty stacks
                continue
            for to_stack in range(len(current_state)):
                if from_stack == to_stack:
                    continue
                # Generate new state
                new_state, block = move_block(current_state, from_stack, to_stack)
                if new_state not in visited:
                    visited.add(new_state)
                    new_moves = moves + [f"Move {block} from {from_stack + 1} to {to_stack + 1}"]
                    queue.append((new_state, new_moves))

    return None

# Initial and goal states
initial_state = (('D', 'I'), ('H', 'A', 'C'), ('F', 'B', 'E'), ('G',))
goal_state = (('A', 'B', 'C', 'D', 'F'), ('E', 'H', 'I'), ('G',))

# Perform BFS to find the solution
solution_moves = bfs(initial_state, goal_state)

# Print the moves
if solution_moves:
    print("<<<" + "\n".join(solution_moves) + ">>>")
else:
    print("No solution found.")