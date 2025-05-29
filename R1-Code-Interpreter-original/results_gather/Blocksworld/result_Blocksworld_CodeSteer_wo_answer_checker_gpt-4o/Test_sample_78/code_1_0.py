from collections import deque
import time

def is_goal_state(state, goal_state):
    return state == goal_state

def get_possible_moves(state):
    moves = []
    for i, stack in enumerate(state):
        if stack:  # If the stack is not empty
            block = stack[-1]
            for j in range(len(state)):
                if i != j:  # Move to a different stack
                    new_state = [list(s) for s in state]  # Deep copy
                    new_state[i].pop()
                    new_state[j].append(block)
                    moves.append((tuple(tuple(s) for s in new_state), f"Move {block} from {i+1} to {j+1}"))
    return moves

def bfs(initial_state, goal_state, timeout=10):
    start_time = time.time()
    queue = deque([(initial_state, [])])
    visited = set()

    while queue:
        if time.time() - start_time > timeout:
            print("Timeout reached. Exiting search.")
            return None

        current_state, path = queue.popleft()
        if current_state in visited:
            continue

        visited.add(current_state)

        if is_goal_state(current_state, goal_state):
            return path

        for new_state, move in get_possible_moves(current_state):
            if new_state not in visited:
                queue.append((new_state, path + [move]))

    return None

initial_state = (('F', 'I', 'H'), ('G', 'A', 'E', 'C', 'B'), ('D',))
goal_state = (('A', 'C'), ('E', 'F', 'G'), ('B', 'D', 'H', 'I'))

solution = bfs(initial_state, goal_state)

if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")