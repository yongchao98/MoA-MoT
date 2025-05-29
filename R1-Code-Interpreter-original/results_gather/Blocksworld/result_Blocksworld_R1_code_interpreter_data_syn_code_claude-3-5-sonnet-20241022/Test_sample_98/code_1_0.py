def find_solution(initial_stacks, goal_stacks):
    def get_top_block(stack):
        return stack[-1] if stack else None
        
    def make_move(stacks, from_stack, to_stack):
        new_stacks = [list(s) for s in stacks]
        block = new_stacks[from_stack].pop()
        new_stacks[to_stack].append(block)
        return new_stacks
        
    def state_to_string(stacks):
        return ';'.join(','.join(stack) for stack in stacks)
        
    def is_valid_move(stacks, from_stack, to_stack):
        if not stacks[from_stack]:  # source stack is empty
            return False
        return True
        
    def is_goal_state(current, goal):
        return all(current[i] == goal[i] for i in range(len(current)))
        
    def find_path(current_stacks, visited, path):
        if is_goal_state(current_stacks, goal_stacks):
            return path
            
        current_state = state_to_string(current_stacks)
        if current_state in visited:
            return None
            
        visited.add(current_state)
        
        # Try all possible moves
        for from_stack in range(len(current_stacks)):
            for to_stack in range(len(current_stacks)):
                if from_stack != to_stack and is_valid_move(current_stacks, from_stack, to_stack):
                    # Get the block to be moved
                    block = current_stacks[from_stack][-1]
                    
                    # Make the move
                    new_stacks = make_move(current_stacks, from_stack, to_stack)
                    
                    # Add move to path
                    new_path = path + [(block, from_stack + 1, to_stack + 1)]
                    
                    # Recursive call
                    result = find_path(new_stacks, visited, new_path)
                    if result is not None:
                        return result
                        
        return None

    # Find solution
    solution = find_path(initial_stacks, set(), [])
    
    if solution:
        # Format and print solution
        moves = []
        for block, source, dest in solution:
            moves.append(f"Move {block} from {source} to {dest}")
        print("<<<" + "\n".join(moves) + ">>>")
    else:
        print("No solution found")

# Initial and goal states
initial = [['A', 'E'], ['D'], ['C', 'B']]
goal = [['A', 'B', 'E'], ['C'], ['D']]

find_solution(initial, goal)