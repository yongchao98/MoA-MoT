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
        if result is not None:
            return result
    return None

initial_state = [['B', 'A'], ['C', 'F', 'G'], ['D', 'I', 'E', 'J', 'H']]
goal_state = [['A', 'E', 'I'], ['B', 'C', 'J'], ['D', 'F', 'G', 'H']]

solution = dfs(initial_state, goal_state, [])

if solution:
    print("<<<")
    for block, src, dest in solution:
        print(f"Move {block} from stack{src+1} to stack{dest+1}")
    print(">>>")
else:
    print("No solution found.")