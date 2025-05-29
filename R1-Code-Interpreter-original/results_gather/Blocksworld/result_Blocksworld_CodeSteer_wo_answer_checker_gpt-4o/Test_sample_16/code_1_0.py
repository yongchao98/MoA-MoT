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
                    moves.append((tuple(map(tuple, new_state)), f"Move {block} from {i+1} to {j+1}"))
    return moves

def bfs(initial_state, goal_state):
    queue = deque([(initial_state, [])])
    visited = set()
    visited.add(initial_state)

    while queue:
        current_state, path = queue.popleft()

        if is_goal(current_state, goal_state):
            return path

        for new_state, move in get_possible_moves(current_state):
            if new_state not in visited:
                visited.add(new_state)
                queue.append((new_state, path + [move]))

    return None

# Initial and goal states
initial_state = (('C', 'A', 'F'), ('B',), ('D', 'E'))
goal_state = (('B',), ('A', 'C', 'D'), ('E', 'F'))

# Find the solution
solution = bfs(initial_state, goal_state)

# Print the solution
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")