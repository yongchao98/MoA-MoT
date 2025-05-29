from collections import deque

def is_goal_state(state, goal_state):
    return state == goal_state

def get_possible_moves(state):
    moves = []
    for i, stack in enumerate(state):
        if stack:  # If the stack is not empty
            block = stack[-1]  # Get the top block
            for j in range(len(state)):
                if i != j:  # Move to a different stack
                    new_state = list(map(list, state))
                    new_state[i].pop()  # Remove block from current stack
                    new_state[j].append(block)  # Add block to new stack
                    moves.append((tuple(map(tuple, new_state)), f"Move {block} from {i+1} to {j+1}"))
    return moves

def bfs(initial_state, goal_state):
    queue = deque([(initial_state, [])])
    visited = set()
    while queue:
        current_state, path = queue.popleft()
        if is_goal_state(current_state, goal_state):
            return path
        visited.add(current_state)
        for new_state, move in get_possible_moves(current_state):
            if new_state not in visited:
                queue.append((new_state, path + [move]))

initial_state = (('D',), ('G', 'C', 'A', 'F'), ('B', 'I', 'E', 'H'))
goal_state = (('A', 'F', 'H'), ('C', 'G', 'I'), ('B', 'D', 'E'))

solution = bfs(initial_state, goal_state)
print("<<<" + "\n".join(solution) + ">>>")