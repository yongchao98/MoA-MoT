from collections import deque

def move_blocks_bfs(initial_state, goal_state):
    queue = deque([(initial_state, [])])
    visited = set()

    while queue:
        current_state, moves = queue.popleft()
        state_tuple = tuple(tuple(stack) for stack in current_state)

        if state_tuple in visited:
            continue
        visited.add(state_tuple)

        if current_state == goal_state:
            return moves

        for i, stack in enumerate(current_state):
            if stack:  # If the stack is not empty
                block = stack[-1]  # Get the top block
                for j in range(len(current_state)):
                    if i != j:  # Move to a different stack
                        # Create a new state by moving the block
                        new_state = [list(s) for s in current_state]
                        new_state[i].pop()
                        new_state[j].append(block)
                        # Add the new state to the queue
                        queue.append((new_state, moves + [f"Move {block} from stack{i+1} to stack{j+1}"]))

    return None

initial_state = [[], ['F'], ['G', 'C'], ['A', 'B', 'D', 'E', 'H']]
goal_state = [['B', 'D', 'E'], ['G', 'H'], ['A', 'C', 'F'], []]

solution = move_blocks_bfs(initial_state, goal_state)
if solution:
    print("<<<")
    for move in solution:
        print(move)
    print(">>>")
else:
    print("No solution found.")