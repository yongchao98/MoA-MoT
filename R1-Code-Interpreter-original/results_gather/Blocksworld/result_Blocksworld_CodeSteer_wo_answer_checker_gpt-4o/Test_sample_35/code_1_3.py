def is_goal_state(stacks, goal):
    return stacks == goal

def get_possible_moves(stacks):
    moves = []
    for i, stack in enumerate(stacks):
        if stack:  # If the stack is not empty
            block = stack[-1]
            for j in range(len(stacks)):
                if i != j:  # Move to a different stack
                    new_stacks = [list(s) for s in stacks]  # Deep copy
                    new_stacks[i].pop()
                    new_stacks[j].append(block)
                    moves.append((block, i + 1, j + 1, new_stacks))
    return moves

def dfs(stacks, goal, path, visited):
    if is_goal_state(stacks, goal):
        return path
    visited.add(tuple(tuple(s) for s in stacks))
    for block, src, dest, new_stacks in get_possible_moves(stacks):
        if tuple(tuple(s) for s in new_stacks) not in visited:
            result = dfs(new_stacks, goal, path + [f"Move {block} from stack{src} to stack{dest}"], visited)
            if result:
                return result
    return None

initial_state = [['H', 'G', 'A', 'C'], ['E'], ['F'], ['B', 'D']]
goal_state = [['A', 'B', 'E', 'G'], ['C', 'F'], ['D', 'H']]

solution = dfs(initial_state, goal_state, [], set())
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")