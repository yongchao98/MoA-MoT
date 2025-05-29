def get_top_block(stack):
    return stack[-1] if stack else None

def make_move(stacks, move_from, move_to):
    if not stacks[move_from]:
        return False
    block = stacks[move_from].pop()
    stacks[move_to].append(block)
    return True

def solve_blocksworld():
    # Initial state
    stacks = [
        ['B', 'K', 'J', 'C'],
        ['F', 'H'],
        ['D', 'E', 'G', 'I', 'A']
    ]
    
    # Goal state
    goal = [
        ['C', 'D', 'F', 'I', 'K'],
        ['B', 'H', 'J'],
        ['A', 'E', 'G']
    ]
    
    moves = []
    
    # Helper function to find a block's location
    def find_block(block):
        for i, stack in enumerate(stacks):
            if block in stack:
                return i
        return None
    
    # Helper function to check if a block is at its goal position
    def is_block_in_position(block):
        current_stack_idx = find_block(block)
        goal_stack_idx = None
        for i, stack in enumerate(goal):
            if block in stack:
                goal_stack_idx = i
                break
        
        if current_stack_idx == goal_stack_idx:
            current_stack = stacks[current_stack_idx]
            goal_stack = goal[goal_stack_idx]
            block_idx = current_stack.index(block)
            return (block_idx == goal_stack.index(block) and 
                    all(current_stack[i] == goal_stack[i] for i in range(block_idx)))
        return False

    # Process each stack in the goal state from bottom to top
    for target_stack_idx, target_stack in enumerate(goal):
        for target_pos, block in enumerate(target_stack):
            if not is_block_in_position(block):
                # Find where the block currently is
                current_stack_idx = find_block(block)
                current_stack = stacks[current_stack_idx]
                
                # Clear blocks above the target block
                while get_top_block(current_stack) != block:
                    top_block = current_stack[-1]
                    # Move to any available stack that's not the current or target
                    temp_stack_idx = next(i for i in range(3) 
                                        if i != current_stack_idx and i != target_stack_idx)
                    make_move(stacks, current_stack_idx, temp_stack_idx)
                    moves.append(f"Move {top_block} from {current_stack_idx + 1} to {temp_stack_idx + 1}")
                
                # Move the block to its target position
                make_move(stacks, current_stack_idx, target_stack_idx)
                moves.append(f"Move {block} from {current_stack_idx + 1} to {target_stack_idx + 1}")

    print('\n'.join(moves))

solve_blocksworld()