def dfs(current_state, goal_state, path):
    if current_state == goal_state:
        return path
    
    for i, stack in enumerate(current_state):
        if stack:  # If the stack is not empty
            block = stack[-1]  # Get the top block
            for j in range(3):
                if i != j:  # Move to a different stack
                    # Create a new state by moving the block
                    new_state = [list(s) for s in current_state]
                    new_state[i].pop()
                    new_state[j].append(block)
                    
                    # Add the move to the path
                    new_path = path + [f"Move {block} from stack{i+1} to stack{j+1}"]
                    
                    # Recursively search from the new state
                    result = dfs(new_state, goal_state, new_path)
                    if result:
                        return result
    return None

initial_state = [['B', 'G'], ['F', 'E', 'A', 'C'], ['H', 'D']]
goal_state = [['A', 'B', 'F'], ['D', 'E'], ['C', 'G', 'H']]

solution = dfs(initial_state, goal_state, [])
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")