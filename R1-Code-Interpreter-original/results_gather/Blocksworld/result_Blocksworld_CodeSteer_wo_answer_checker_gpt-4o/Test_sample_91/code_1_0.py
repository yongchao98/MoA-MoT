from collections import deque

def is_goal(state):
    # Check if the current state matches the goal state
    return state == (('D', 'E', 'G', 'H'), ('B', 'C', 'F'), ('A',), ())

def get_possible_moves(state):
    moves = []
    for i, stack in enumerate(state):
        if stack:  # If the stack is not empty
            block = stack[-1]  # Get the top block
            for j in range(len(state)):
                if i != j:  # Move to a different stack
                    new_state = [list(s) for s in state]
                    new_state[i].pop()
                    new_state[j].append(block)
                    moves.append((f"Move {block} from {i+1} to {j+1}", tuple(map(tuple, new_state))))
    return moves

def bfs(initial_state):
    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(initial_state)

    while queue:
        current_state, path = queue.popleft()

        # Debugging: Print the current state and path
        print(f"Current state: {current_state}, Path: {path}")

        if is_goal(current_state):
            return path

        for move, new_state in get_possible_moves(current_state):
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + [move]))

initial_state = (('E',), ('B',), ('H', 'C', 'A', 'F'), ('D', 'G'))
solution = bfs(initial_state)

if solution:
    print("<<<")
    for move in solution:
        print(move)
    print(">>>")
else:
    print("No solution found.")