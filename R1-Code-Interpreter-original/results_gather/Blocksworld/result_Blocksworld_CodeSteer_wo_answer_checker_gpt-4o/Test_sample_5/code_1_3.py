from collections import deque

def is_goal(state, goal):
    return state == goal

def get_possible_moves(state):
    moves = []
    for i, stack in enumerate(state):
        if stack:  # If the stack is not empty
            block = stack[-1]
            for j in range(len(state)):
                if i != j:  # Move to a different stack
                    new_state = [s[:] for s in state]  # Deep copy of the state
                    new_state[i].pop()
                    new_state[j].append(block)
                    moves.append((new_state, f"Move {block} from {i+1} to {j+1}"))
    return moves

def bfs(initial_state, goal_state):
    queue = deque([(initial_state, [])])
    visited = set()

    while queue:
        current_state, path = queue.popleft()
        state_tuple = tuple(tuple(stack) for stack in current_state)

        if state_tuple in visited:
            continue

        visited.add(state_tuple)

        if is_goal(current_state, goal_state):
            return path

        for new_state, move in get_possible_moves(current_state):
            queue.append((new_state, path + [move]))

    return None

initial_state = [['D', 'C', 'F', 'B'], ['E'], [], ['G', 'A']]
goal_state = [['B', 'F', 'G'], ['C', 'E'], ['A', 'D'], []]

solution = bfs(initial_state, goal_state)

if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")