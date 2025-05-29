def print_moves(initial_stacks, goal_stacks):
    moves = []
    current_stacks = [list(stack) for stack in initial_stacks]
    
    # Helper function to find a block's location
    def find_block(block):
        for i, stack in enumerate(current_stacks):
            if block in stack:
                return i
        return None
    
    # Helper function to get top block
    def get_top(stack_idx):
        return current_stacks[stack_idx][-1] if current_stacks[stack_idx] else None
    
    # Helper function to make a move
    def make_move(block, from_stack, to_stack):
        if get_top(from_stack) == block:
            current_stacks[from_stack].pop()
            current_stacks[to_stack].append(block)
            moves.append(f"Move {block} from {from_stack + 1} to {to_stack + 1}")
            return True
        return False

    # Process each stack in goal state
    for goal_stack_idx, goal_stack in enumerate(goal_stacks):
        for block in goal_stack:
            current_loc = find_block(block)
            if current_loc is None:
                continue
                
            # If block is not in correct position
            if current_loc != goal_stack_idx or current_stacks[current_loc].index(block) != len(current_stacks[current_loc]) - 1:
                # Move blocks above it to other stacks temporarily
                while get_top(current_loc) != block:
                    top_block = get_top(current_loc)
                    # Find a temporary stack
                    for temp_stack in range(3):
                        if temp_stack != current_loc and temp_stack != goal_stack_idx:
                            make_move(top_block, current_loc, temp_stack)
                            break
                
                # Move the block to its goal position
                make_move(block, current_loc, goal_stack_idx)

    print('<<<')
    print('\n'.join(moves))
    print('>>>')

# Initial and goal states
initial = [['G', 'C'], ['D', 'A', 'E', 'B'], ['I', 'H', 'F']]
goal = [['B', 'I'], ['E'], ['A', 'C', 'D', 'F', 'G', 'H']]

print_moves(initial, goal)