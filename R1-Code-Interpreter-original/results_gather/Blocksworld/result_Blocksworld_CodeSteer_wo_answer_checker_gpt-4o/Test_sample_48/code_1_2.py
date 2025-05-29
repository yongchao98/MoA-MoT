def dfs(current_state, goal_state, path, visited, depth, max_depth):
    # Base case: check if the current state matches the goal state
    if current_state == goal_state:
        return path
    
    # Check if the current state has been visited or if the depth limit is reached
    state_tuple = tuple(tuple(stack) for stack in current_state)
    if state_tuple in visited or depth > max_depth:
        return None
    
    # Mark the current state as visited
    visited.add(state_tuple)
    
    # Debugging: print the current state and path
    print(f"Depth: {depth}, Current State: {current_state}, Path: {path}")
    
    # Explore possible moves
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
                    result = dfs(new_state, goal_state, new_path, visited, depth + 1, max_depth)
                    if result:
                        return result
    
    # Unmark the current state as visited (backtrack)
    visited.remove(state_tuple)
    return None

initial_state = [['B', 'G'], ['F', 'E', 'A', 'C'], ['H', 'D']]
goal_state = [['A', 'B', 'F'], ['D', 'E'], ['C', 'G', 'H']]

# Set a reasonable depth limit to prevent recursion errors
max_depth = 20
solution = dfs(initial_state, goal_state, [], set(), 0, max_depth)
if solution:
    print("<<<" + "\n".join(solution) + ">>>")
else:
    print("No solution found.")