from collections import deque

def bfs(initial_state, goal_state):
    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(tuple(tuple(stack) for stack in initial_state))

    while queue:
        current_state, path = queue.popleft()

        if current_state == goal_state:
            return path

        for i, stack in enumerate(current_state):
            if stack:  # If the stack is not empty
                block = stack[-1]  # Get the top block
                for j in range(len(current_state)):
                    if i != j:  # Move to a different stack
                        # Create a new state by moving the block
                        new_state = [list(s) for s in current_state]
                        new_state[i].pop()
                        new_state[j].append(block)
                        new_state_tuple = tuple(tuple(stack) for stack in new_state)
                        if new_state_tuple not in visited:
                            visited.add(new_state_tuple)
                            move = f"Move {block} from stack{i+1} to stack{j+1}"
                            queue.append((new_state, path + [move]))

    return None

initial_state = [['C', 'E'], ['D'], ['G', 'F', 'B'], ['A', 'H']]
goal_state = [['B', 'C', 'D', 'E', 'F', 'G', 'H'], [], ['A'], []]

solution = bfs(initial_state, goal_state)

if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")