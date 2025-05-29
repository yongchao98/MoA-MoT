def is_goal_state(stacks, goal):
    return stacks == goal

def dfs(stacks, goal, path, visited):
    if is_goal_state(stacks, goal):
        return path

    visited.add(tuple(tuple(stack) for stack in stacks))

    for i in range(len(stacks)):
        if not stacks[i]:
            continue
        block = stacks[i][-1]
        for j in range(len(stacks)):
            if i != j:
                new_stacks = [stack[:] for stack in stacks]
                new_stacks[i].pop()
                new_stacks[j].append(block)
                state_tuple = tuple(tuple(stack) for stack in new_stacks)
                if state_tuple not in visited:
                    result = dfs(new_stacks, goal, path + [f"Move {block} from stack{i+1} to stack{j+1}"], visited)
                    if result:
                        return result
    return None

initial_state = [['B', 'H', 'G', 'I', 'D'], ['A', 'F'], ['E'], ['J', 'C']]
goal_state = [['H'], ['A', 'B', 'D', 'G'], ['C', 'E', 'F', 'I', 'J'], []]

solution = dfs(initial_state, goal_state, [], set())
if solution:
    print("<<<")
    for move in solution:
        print(move)
    print(">>>")
else:
    print("No solution found.")