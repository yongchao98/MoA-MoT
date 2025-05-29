def get_top_blocks(state):
    return [stack[-1] if stack else None for stack in state]

def is_valid_move(state, from_stack, to_stack, goal_state):
    if not state[from_stack]:  # Source stack is empty
        return False
    
    moving_block = state[from_stack][-1]
    
    # Check if this move gets us closer to the goal
    if to_stack < len(goal_state):
        if not goal_state[to_stack]:  # If goal stack is empty
            return True
        if not state[to_stack]:  # If destination stack is empty
            return moving_block in goal_state[to_stack]
        # Check if this block should be in this stack and should be above the current top block
        return (moving_block in goal_state[to_stack] and 
                state[to_stack][-1] in goal_state[to_stack] and 
                goal_state[to_stack].index(moving_block) > goal_state[to_stack].index(state[to_stack][-1]))
    return True

def solve_blocksworld():
    # Initial and goal states
    initial_state = [['J', 'B'], ['F', 'I', 'E'], ['D', 'A', 'C', 'H', 'G']]
    goal_state = [['B', 'C', 'G', 'H', 'I', 'J'], ['D', 'F'], ['A', 'E']]
    
    moves = []
    current_state = [list(stack) for stack in initial_state]
    
    # Helper function to find block location
    def find_block(block):
        for i, stack in enumerate(current_state):
            if block in stack:
                return i
        return None
    
    # Helper function to clear blocks above target
    def clear_top_block(stack_idx):
        while current_state[stack_idx]:
            top_block = current_state[stack_idx][-1]
            # Find where this block should go
            for goal_idx, goal_stack in enumerate(goal_state):
                if top_block in goal_stack:
                    if is_valid_move(current_state, stack_idx, goal_idx, goal_state):
                        current_state[goal_idx].append(current_state[stack_idx].pop())
                        moves.append(f"Move {top_block} from {stack_idx + 1} to {goal_idx + 1}")
                        return True
            # If can't move to goal position yet, move to temporary position
            for temp_idx in range(3):
                if temp_idx != stack_idx and is_valid_move(current_state, stack_idx, temp_idx, goal_state):
                    current_state[temp_idx].append(current_state[stack_idx].pop())
                    moves.append(f"Move {top_block} from {stack_idx + 1} to {temp_idx + 1}")
                    return True
        return False

    # Main solving loop
    for target_stack_idx, target_stack in enumerate(goal_state):
        for target_pos, target_block in enumerate(target_stack):
            current_loc = find_block(target_block)
            if current_loc is None:
                continue
                
            # If block is not in correct position
            if (current_loc != target_stack_idx or 
                current_state[current_loc].index(target_block) != target_pos):
                
                # Clear blocks above target block
                while current_state[current_loc][-1] != target_block:
                    clear_top_block(current_loc)
                
                # Move block to target position
                if is_valid_move(current_state, current_loc, target_stack_idx, goal_state):
                    current_state[target_stack_idx].append(current_state[current_loc].pop())
                    moves.append(f"Move {target_block} from {current_loc + 1} to {target_stack_idx + 1}")

    return moves

# Get and print solution
solution = solve_blocksworld()
print("<<<" + "\n".join(solution) + ">>>")