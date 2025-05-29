from collections import deque

def bfs(initial_state, goal_state):
    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(tuple(map(tuple, initial_state)))

    while queue:
        current_state, path = queue.popleft()

        if current_state == goal_state:
            return path

        for i in range(3):
            if current_state[i]:  # If the stack is not empty
                block = current_state[i][-1]
                for j in range(3):
                    if i != j:  # Move to a different stack
                        new_state = [list(stack) for stack in current_state]
                        new_state[i].pop()
                        new_state[j].append(block)

                        state_tuple = tuple(map(tuple, new_state))
                        if state_tuple not in visited:
                            visited.add(state_tuple)
                            queue.append((new_state, path + [f"Move {block} from {i+1} to {j+1}"]))

    return None

initial_state = [['F', 'B', 'I', 'K', 'A'], ['H', 'D'], ['J', 'E', 'C', 'G']]
goal_state = [['C', 'G'], ['A', 'B', 'D', 'E', 'I', 'K'], ['F', 'H', 'J']]

solution = bfs(initial_state, goal_state)
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")