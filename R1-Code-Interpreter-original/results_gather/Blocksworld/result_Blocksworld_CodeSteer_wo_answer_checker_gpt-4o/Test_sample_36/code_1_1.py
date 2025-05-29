def dfs(initial_state, goal_state, path, visited):
    if initial_state == goal_state:
        return path

    visited.add(tuple(map(tuple, initial_state)))

    for i in range(len(initial_state)):
        if not initial_state[i]:
            continue
        block = initial_state[i][-1]
        for j in range(len(initial_state)):
            if i != j:
                # Make a move
                new_state = [stack[:] for stack in initial_state]
                new_state[i].pop()
                new_state[j].append(block)
                state_tuple = tuple(map(tuple, new_state))
                if state_tuple not in visited:
                    result = dfs(new_state, goal_state, path + [f"Move {block} from stack{i+1} to stack{j+1}"], visited)
                    if result:
                        return result
    return None

initial_state = [['D', 'E'], ['G', 'I', 'F'], ['C', 'B', 'H', 'A']]
goal_state = [['B', 'C', 'E', 'F', 'I'], ['G', 'H'], ['A', 'D']]
visited = set()
solution = dfs(initial_state, goal_state, [], visited)

if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")