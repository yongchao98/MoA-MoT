def get_top_block(stack):
    return stack[-1] if stack else None

def is_valid_move(stacks, from_stack, to_stack):
    return bool(stacks[from_stack-1])

def make_move(stacks, from_stack, to_stack):
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    return f"Move {block} from {from_stack} to {to_stack}"

def find_block_location(stacks, block):
    for i, stack in enumerate(stacks):
        if stack and block in stack:
            return i + 1
    return None

def solve_blocksworld(current, goal):
    moves = []
    current = [list(stack) for stack in current]
    
    # Process each goal stack from bottom to top
    for goal_stack_idx, goal_stack in enumerate(goal):
        for height, target_block in enumerate(goal_stack):
            # Skip if block is already in correct position
            if (len(current[goal_stack_idx]) > height and 
                current[goal_stack_idx][height] == target_block):
                continue
                
            # Find where the target block currently is
            source_stack = find_block_location(current, target_block)
            if source_stack is None:
                continue
                
            # If target block is not on top, move blocks above it to other stacks
            while get_top_block(current[source_stack-1]) != target_block:
                top_block = get_top_block(current[source_stack-1])
                # Move to any available stack except source and destination
                for temp_stack in range(1, 4):
                    if temp_stack != source_stack and temp_stack != goal_stack_idx + 1:
                        moves.append(make_move(current, source_stack, temp_stack))
                        break
            
            # Now move the target block to its goal position
            if get_top_block(current[source_stack-1]) == target_block:
                # If destination stack has blocks above target position, move them temporarily
                while len(current[goal_stack_idx]) > height:
                    top_block = current[goal_stack_idx].pop()
                    for temp_stack in range(1, 4):
                        if temp_stack != source_stack and temp_stack != goal_stack_idx + 1:
                            current[temp_stack-1].append(top_block)
                            moves.append(f"Move {top_block} from {goal_stack_idx + 1} to {temp_stack}")
                            break
                
                # Move target block to goal position
                moves.append(make_move(current, source_stack, goal_stack_idx + 1))

    return moves

# Initial and goal states
initial = [['F'], ['G', 'B', 'D'], ['A', 'I', 'C', 'E', 'H']]
goal = [['E', 'G', 'H'], ['A', 'D', 'I'], ['B', 'C', 'F']]

# Find and print solution
solution = solve_blocksworld(initial, goal)
print('<<<')
print('\n'.join(solution))
print('>>>')