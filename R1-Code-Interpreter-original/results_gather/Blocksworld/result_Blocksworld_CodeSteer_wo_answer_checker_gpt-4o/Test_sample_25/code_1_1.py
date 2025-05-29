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

def dfs(stacks, goal, path):
    if is_goal_state(stacks, goal):
        return path
    for move in get_possible_moves(stacks):
        new_stacks = apply_move(stacks, move)
        result = dfs(new_stacks, goal, path + [move])
        if result:
            return result
    return None

initial_state = [['A', 'G'], ['F', 'E'], ['C'], ['H', 'J', 'I', 'B', 'D']]
goal_state = [['A', 'B'], ['C', 'E', 'G', 'J'], ['D', 'F', 'H', 'I']]

solution = dfs(initial_state, goal_state, [])
if solution:
    print("<<<" + "\n".join(f"Move {block} from {src+1} to {dest+1}" for block, src, dest in solution) + ">>>")
else:
    print("No solution found.")