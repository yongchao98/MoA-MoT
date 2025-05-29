from collections import deque

def is_goal(state, goal):
    return state == goal

def get_possible_moves(state):
    moves = []
    for i, stack in enumerate(state):
        if stack:  # If the stack is not empty
            block = stack[-1]  # Get the top block
            for j in range(len(state)):
                if i != j:  # Move to a different stack
                    new_state = list(map(list, state))
                    new_state[i].pop()
                    new_state[j].append(block)
                    moves.append((block, i + 1, j + 1, tuple(map(tuple, new_state))))
    return moves

def bfs(initial_state, goal_state):
    queue = deque([(initial_state, [])])
    visited = set()

    while queue:
        current_state, path = queue.popleft()
        if is_goal(current_state, goal_state):
            return path

        if current_state in visited:
            continue

        visited.add(current_state)

        for block, src, dest, new_state in get_possible_moves(current_state):
            if new_state not in visited:
                queue.append((new_state, path + [f"Move {block} from {src} to {dest}"]))

    return None

initial_state = (('F',), ('B',), ('D', 'A', 'E', 'C'))
goal_state = (('F',), ('C', 'E'), ('A', 'B', 'D'))

solution = bfs(initial_state, goal_state)

if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")