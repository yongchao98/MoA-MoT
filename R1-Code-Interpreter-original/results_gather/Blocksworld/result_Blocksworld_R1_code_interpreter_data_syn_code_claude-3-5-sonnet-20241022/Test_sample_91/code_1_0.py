def get_state_string(stacks):
    return ';'.join(','.join(stack) for stack in stacks)

def get_moves():
    # Initial and goal states
    initial = [['E'], ['B'], ['H', 'C', 'A', 'F'], ['D', 'G']]
    goal = [['D', 'E', 'G', 'H'], ['B', 'C', 'F'], ['A'], []]
    goal_str = get_state_string(goal)
    
    # Stack for DFS [(current_state, moves_so_far)]
    stack = [(initial, [])]
    visited = set()
    max_depth = 12  # Limit search depth
    
    while stack:
        current, moves = stack.pop()
        
        # Skip if too many moves
        if len(moves) > max_depth:
            continue
            
        current_str = get_state_string(current)
        if current_str in visited:
            continue
            
        visited.add(current_str)
        
        # Check if goal reached
        if current_str == goal_str:
            return moves
            
        # Generate possible moves
        for i in range(len(current)):
            if not current[i]:  # Skip empty stacks
                continue
            block = current[i][-1]  # Get top block
            
            # Try moving to each other stack
            for j in range(len(current)):
                if i != j:
                    # Make move
                    new_state = [stack[:] for stack in current]
                    new_state[i].pop()
                    new_state[j].append(block)
                    
                    move = f"Move {block} from {i+1} to {j+1}"
                    stack.append((new_state, moves + [move]))

# Run solution and print result
solution = get_moves()
print("<<<" + "\n".join(solution) + ">>>")