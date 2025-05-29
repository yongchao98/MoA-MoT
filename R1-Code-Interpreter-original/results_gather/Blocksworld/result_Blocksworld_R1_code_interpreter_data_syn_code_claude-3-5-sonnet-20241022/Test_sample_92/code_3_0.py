def get_top_block(stack):
    return stack[-1] if stack else None

def make_move(stacks, from_stack, to_stack, moves):
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    moves.append(f"Move {block} from {from_stack} to {to_stack}")

def find_block_location(stacks, block):
    for i, stack in enumerate(stacks):
        if block in stack:
            return i + 1, stack.index(block)
    return None, None

def can_move_block(stacks, block):
    stack_num, pos = find_block_location(stacks, block)
    if stack_num is None:
        return False
    return pos == len(stacks[stack_num-1]) - 1

def solve_blocksworld(initial, goal):
    current = [list(stack) for stack in initial]
    moves = []
    
    # Build goal stacks one by one
    for goal_stack_idx, goal_stack in enumerate(goal):
        for target_pos, target_block in enumerate(goal_stack):
            # Skip if block is already in correct position
            if (len(current[goal_stack_idx]) > target_pos and 
                current[goal_stack_idx][target_pos] == target_block):
                continue
            
            # Clear the destination stack if needed
            while len(current[goal_stack_idx]) > target_pos:
                temp_block = current[goal_stack_idx][-1]
                # Find a temporary stack
                for temp_stack in range(3):
                    if temp_stack != goal_stack_idx and temp_stack != find_block_location(current, target_block)[0] - 1:
                        make_move(current, goal_stack_idx + 1, temp_stack + 1, moves)
                        break
            
            # Keep clearing blocks until target block is accessible
            while not can_move_block(current, target_block):
                src_stack, _ = find_block_location(current, target_block)
                top_block = get_top_block(current[src_stack-1])
                # Move blocking block to any available stack
                for temp_stack in range(3):
                    if temp_stack != src_stack - 1 and temp_stack != goal_stack_idx:
                        make_move(current, src_stack, temp_stack + 1, moves)
                        break
            
            # Move target block to goal position
            src_stack, _ = find_block_location(current, target_block)
            make_move(current, src_stack, goal_stack_idx + 1, moves)

    return moves

# Initial and goal states
initial = [['F'], ['G', 'B', 'D'], ['A', 'I', 'C', 'E', 'H']]
goal = [['E', 'G', 'H'], ['A', 'D', 'I'], ['B', 'C', 'F']]

# Find and print solution
solution = solve_blocksworld(initial, goal)
print('<<<')
print('\n'.join(solution))
print('>>>')