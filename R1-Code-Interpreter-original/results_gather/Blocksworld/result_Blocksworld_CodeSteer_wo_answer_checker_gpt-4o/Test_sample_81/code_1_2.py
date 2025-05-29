from collections import deque

def is_goal_state(state, goal_state):
    return state == goal_state

def get_valid_moves(state):
    moves = []
    for i, stack in enumerate(state):
        if stack:  # If the stack is not empty
            block = stack[-1]
            for j in range(len(state)):
                if i != j:  # Move to a different stack
                    moves.append((block, i, j))
    return moves

def apply_move(state, move):
    block, src, dest = move
    new_state = [stack[:] for stack in state]  # Deep copy of the state
    new_state[src].pop()
    new_state[dest].append(block)
    return new_state

def bfs(initial_state, goal_state):
    queue = deque([(initial_state, [])])
    visited = set()

    while queue:
        current_state, path = queue.popleft()

        if is_goal_state(current_state, goal_state):
            return path

        state_tuple = tuple(tuple(stack) for stack in current_state)
        if state_tuple in visited:
            continue

        visited.add(state_tuple)

        for move in get_valid_moves(current_state):
            new_state = apply_move(current_state, move)
            queue.append((new_state, path + [move]))

    return None

initial_state = [['D'], ['A', 'J'], ['H', 'F', 'C', 'E', 'I'], ['B', 'G']]
goal_state = [['F', 'I'], ['A', 'B', 'C', 'D', 'G'], ['E', 'H', 'J']]

solution = bfs(initial_state, goal_state)

if solution:
    print("<<<" + "\n".join(f"Move {block} from stack{src+1} to stack{dest+1}" for block, src, dest in solution) + ">>>")
else:
    print("No solution found.")