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

def dfs(state, goal_state, path):
    if is_goal_state(state, goal_state):
        return path

    for move in get_valid_moves(state):
        new_state = apply_move(state, move)
        result = dfs(new_state, goal_state, path + [move])
        if result is not None:
            return result

    return None

initial_state = [['D'], ['A', 'J'], ['H', 'F', 'C', 'E', 'I'], ['B', 'G']]
goal_state = [['F', 'I'], ['A', 'B', 'C', 'D', 'G'], ['E', 'H', 'J']]

solution = dfs(initial_state, goal_state, [])

if solution:
    print("<<<" + "\n".join(f"Move {block} from stack{src+1} to stack{dest+1}" for block, src, dest in solution) + ">>>")
else:
    print("No solution found.")