from collections import deque

def generate_moves(state):
    moves = []
    for i, stack in enumerate(state):
        if stack:  # If the stack is not empty
            block = stack[-1]  # Get the top block
            for j in range(len(state)):
                if i != j:  # Move to a different stack
                    new_state = list(map(list, state))
                    new_state[i].pop()  # Remove block from current stack
                    new_state[j].append(block)  # Add block to new stack
                    moves.append((tuple(map(tuple, new_state)), f"Move {block} from stack{i+1} to stack{j+1}"))
    return moves

def bfs(initial_state, goal_state):
    queue = deque([(initial_state, [])])  # Each element is (state, path)
    visited = set()

    while queue:
        current_state, path = queue.popleft()
        if current_state == goal_state:
            return path

        if current_state in visited:
            continue

        visited.add(current_state)

        for new_state, move in generate_moves(current_state):
            if new_state not in visited:
                queue.append((new_state, path + [move]))

    return None

# Test with the given problem
initial_state = (('C', 'F'), ('D',), ('A', 'B', 'E'))
goal_state = (('E', 'F'), ('A', 'B', 'C', 'D'), ())

solution = bfs(initial_state, goal_state)

if solution:
    print("<<<")
    for move in solution:
        print(move)
    print(">>>")
else:
    print("No solution found.")