def is_goal_state(stacks, goal_stacks):
    return stacks == goal_stacks

def get_possible_moves(stacks):
    moves = []
    for i, stack in enumerate(stacks):
        if stack:  # If the stack is not empty
            block = stack[-1]
            for j, target_stack in enumerate(stacks):
                if i != j:  # Don't move to the same stack
                    moves.append((block, i, j))
    return moves

def apply_move(stacks, move):
    block, source, destination = move
    new_stacks = [stack[:] for stack in stacks]
    new_stacks[source].pop()
    new_stacks[destination].append(block)
    return new_stacks

def dfs_limited(stacks, goal_stacks, path, visited, depth_limit):
    if is_goal_state(stacks, goal_stacks):
        return path
    if depth_limit == 0:
        return None

    for move in get_possible_moves(stacks):
        new_stacks = apply_move(stacks, move)
        stacks_tuple = tuple(tuple(stack) for stack in new_stacks)
        if stacks_tuple not in visited:
            visited.add(stacks_tuple)
            result = dfs_limited(new_stacks, goal_stacks, path + [move], visited, depth_limit - 1)
            if result:
                return result
            visited.remove(stacks_tuple)
    return None

def iddfs(initial_stacks, goal_stacks, max_depth):
    for depth in range(max_depth):
        visited = set()
        initial_tuple = tuple(tuple(stack) for stack in initial_stacks)
        visited.add(initial_tuple)
        result = dfs_limited(initial_stacks, goal_stacks, [], visited, depth)
        if result:
            return result
    return None

initial_stacks = [
    ['I'],
    ['G', 'C', 'H', 'A', 'E'],
    ['F', 'B'],
    ['D', 'K', 'J']
]

goal_stacks = [
    ['A', 'D', 'F', 'I', 'J'],
    ['G', 'H'],
    ['B', 'C', 'E', 'K']
]

max_depth = 20  # Set a reasonable depth limit
solution = iddfs(initial_stacks, goal_stacks, max_depth)

if solution:
    print("<<<")
    for move in solution:
        block, source, destination = move
        print(f"Move {block} from {source + 1} to {destination + 1}")
    print(">>>")
else:
    print("No solution found.")