def is_goal_state(state, goal):
    return state == goal

def get_valid_moves(state):
    moves = []
    for i, stack in enumerate(state):
        if stack:  # If the stack is not empty
            block = stack[-1]  # Get the top block
            for j in range(len(state)):
                if i != j:  # Move to a different stack
                    moves.append((block, i, j))
    return moves

def apply_move(state, move):
    block, src, dest = move
    new_state = [list(stack) for stack in state]  # Deep copy the state
    new_state[src].pop()  # Remove block from source stack
    new_state[dest].append(block)  # Add block to destination stack
    return new_state

def dfs(initial_state, goal_state):
    stack = [(initial_state, [])]  # Stack of (state, moves)
    visited = set()

    while stack:
        current_state, moves = stack.pop()
        state_tuple = tuple(tuple(stack) for stack in current_state)
        
        if state_tuple in visited:
            continue
        visited.add(state_tuple)

        if is_goal_state(current_state, goal_state):
            return moves

        for move in get_valid_moves(current_state):
            new_state = apply_move(current_state, move)
            stack.append((new_state, moves + [move]))

    return None

initial_state = [['C'], ['B', 'E', 'A', 'J', 'F'], ['G', 'I'], ['D', 'H']]
goal_state = [['E'], ['A', 'B', 'C', 'D', 'H'], ['F', 'G', 'I', 'J']]

solution_moves = dfs(initial_state, goal_state)

if solution_moves:
    print("<<<")
    for block, src, dest in solution_moves:
        print(f"Move {block} from {src + 1} to {dest + 1}")
    print(">>>")
else:
    print("No solution found.")