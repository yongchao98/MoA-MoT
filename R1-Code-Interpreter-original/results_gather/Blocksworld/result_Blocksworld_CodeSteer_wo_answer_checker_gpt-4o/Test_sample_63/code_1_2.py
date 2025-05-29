def dfs(current_state, goal_state, path, visited, depth, max_depth):
    if depth > max_depth:
        return None
    if current_state == goal_state:
        return path

    visited.add(tuple(tuple(stack) for stack in current_state))

    for i, stack in enumerate(current_state):
        if stack:  # If the stack is not empty
            block = stack[-1]  # Get the top block
            for j in range(len(current_state)):
                if i != j:  # Move to a different stack
                    # Create a new state by moving the block
                    new_state = [list(s) for s in current_state]
                    new_state[i].pop()
                    new_state[j].append(block)
                    new_state_tuple = tuple(tuple(stack) for stack in new_state)
                    if new_state_tuple not in visited:
                        move = f"Move {block} from stack{i+1} to stack{j+1}"
                        result = dfs(new_state, goal_state, path + [move], visited, depth + 1, max_depth)
                        if result:
                            return result
    return None

initial_state = [['C', 'E'], ['D'], ['G', 'F', 'B'], ['A', 'H']]
goal_state = [['B', 'C', 'D', 'E', 'F', 'G', 'H'], [], ['A'], []]

max_depth = 50  # Set a reasonable maximum depth
visited = set()
solution = dfs(initial_state, goal_state, [], visited, 0, max_depth)

if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")