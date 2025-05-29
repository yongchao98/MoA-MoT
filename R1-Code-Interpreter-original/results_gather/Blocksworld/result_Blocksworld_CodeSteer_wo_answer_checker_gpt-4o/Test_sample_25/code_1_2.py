from collections import deque

def is_goal_state(stacks, goal):
    return stacks == goal

def get_possible_moves(stacks):
    moves = []
    for i, stack in enumerate(stacks):
        if stack:  # If the stack is not empty
            block = stack[-1]
            for j in range(len(stacks)):
                if i != j:  # Move to a different stack
                    moves.append((block, i, j))
    return moves

def apply_move(stacks, move):
    block, src, dest = move
    new_stacks = [list(stack) for stack in stacks]  # Deep copy
    new_stacks[src].pop()
    new_stacks[dest].append(block)
    return new_stacks

def bfs(initial_state, goal_state):
    queue = deque([(initial_state, [])])
    visited = set()

    while queue:
        current_state, path = queue.popleft()
        state_tuple = tuple(tuple(stack) for stack in current_state)

        if state_tuple in visited:
            continue

        visited.add(state_tuple)

        if is_goal_state(current_state, goal_state):
            return path

        for move in get_possible_moves(current_state):
            new_state = apply_move(current_state, move)
            queue.append((new_state, path + [move]))

    return None

initial_state = [['A', 'G'], ['F', 'E'], ['C'], ['H', 'J', 'I', 'B', 'D']]
goal_state = [['A', 'B'], ['C', 'E', 'G', 'J'], ['D', 'F', 'H', 'I']]

solution = bfs(initial_state, goal_state)
if solution:
    print("<<<" + "\n".join(f"Move {block} from {src+1} to {dest+1}" for block, src, dest in solution) + ">>>")
else:
    print("No solution found.")