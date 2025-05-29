def dfs(current_state, goal_state, path):
    if current_state == goal_state:
        return path

    for i, stack in enumerate(current_state):
        if stack:  # If the stack is not empty
            block = stack[-1]  # Get the top block
            for j in range(len(current_state)):
                if i != j:  # Move to a different stack
                    # Create a new state by moving the block
                    new_state = [s[:] for s in current_state]
                    new_state[i].pop()
                    new_state[j].append(block)
                    move = f"Move {block} from stack{i+1} to stack{j+1}"
                    result = dfs(new_state, goal_state, path + [move])
                    if result:
                        return result
    return None

initial_state = [['C', 'E'], ['D'], ['G', 'F', 'B'], ['A', 'H']]
goal_state = [['B', 'C', 'D', 'E', 'F', 'G', 'H'], [], ['A'], []]

solution = dfs(initial_state, goal_state, [])
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")