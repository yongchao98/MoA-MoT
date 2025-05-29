def is_goal_state(stacks, goal):
    return stacks == goal

def get_valid_moves(stacks):
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
    new_stacks = [stack[:] for stack in stacks]
    new_stacks[src].pop()
    new_stacks[dest].append(block)
    return new_stacks

def dfs(stacks, goal, path, visited):
    if is_goal_state(stacks, goal):
        return path
    visited.add(tuple(tuple(stack) for stack in stacks))
    for move in get_valid_moves(stacks):
        new_stacks = apply_move(stacks, move)
        new_state = tuple(tuple(stack) for stack in new_stacks)
        if new_state not in visited:
            result = dfs(new_stacks, goal, path + [move], visited)
            if result is not None:
                return result
    return None

initial_state = [['G', 'A', 'H'], ['E'], ['I', 'C', 'B'], ['F', 'D']]
goal_state = [['A', 'B', 'C', 'D', 'G', 'H', 'I'], ['E'], ['F']]

solution = dfs(initial_state, goal_state, [], set())

if solution:
    print("<<<")
    for block, src, dest in solution:
        print(f"Move {block} from {src + 1} to {dest + 1}")
    print(">>>")
else:
    print("No solution found.")