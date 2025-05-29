def find_solution():
    initial = [['G', 'I', 'D'], ['E', 'A', 'K'], ['C', 'B', 'H', 'J', 'F']]
    goal = [['A', 'D', 'F'], ['B', 'H', 'I', 'J', 'K'], ['C', 'E', 'G']]
    
    def find_block_location(block, stacks):
        for i, stack in enumerate(stacks):
            if stack and block in stack:
                return i, stack.index(block)
        return None
    
    def can_move_block(block, stacks):
        loc = find_block_location(block, stacks)
        if loc:
            stack_idx, block_idx = loc
            return block_idx == len(stacks[stack_idx]) - 1
        return False
    
    current = [list(stack) for stack in initial]
    moves = []
    
    # Helper function to move a block
    def move_block(block, from_stack, to_stack):
        if current[from_stack][-1] == block:
            moves.append(f"Move {block} from {from_stack + 1} to {to_stack + 1}")
            current[from_stack].pop()
            current[to_stack].append(block)
            return True
        return False
    
    # Process each target stack
    for target_stack_idx, target_stack in enumerate(goal):
        for target_pos, target_block in enumerate(target_stack):
            # Find current location of the target block
            curr_loc = find_block_location(target_block, current)
            if curr_loc:
                curr_stack_idx, curr_block_idx = curr_loc
                
                # Clear blocks above the target block
                while curr_block_idx < len(current[curr_stack_idx]) - 1:
                    top_block = current[curr_stack_idx][-1]
                    # Find temporary stack (not current and not target)
                    temp_stack = 3 - curr_stack_idx - target_stack_idx
                    move_block(top_block, curr_stack_idx, temp_stack)
                
                # Move the target block to its destination
                if curr_stack_idx != target_stack_idx:
                    move_block(target_block, curr_stack_idx, target_stack_idx)

    return moves

# Execute and print solution
solution = find_solution()
print('\n'.join(solution))