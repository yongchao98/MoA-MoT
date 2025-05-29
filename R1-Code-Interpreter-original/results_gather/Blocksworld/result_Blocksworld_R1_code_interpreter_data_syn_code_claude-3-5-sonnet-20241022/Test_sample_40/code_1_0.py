def get_top_blocks(stacks):
    return [stack[-1] if stack else None for stack in stacks]

def is_goal_achieved(current, goal):
    return all(len(c) == len(g) and all(x == y for x, y in zip(c, g)) 
              for c, g in zip(current, goal))

def find_solution():
    initial = [['G', 'I', 'D'], ['E', 'A', 'K'], ['C', 'B', 'H', 'J', 'F']]
    goal = [['A', 'D', 'F'], ['B', 'H', 'I', 'J', 'K'], ['C', 'E', 'G']]
    
    moves = []
    current = [list(stack) for stack in initial]
    
    # First clear F, K, J to make space
    while current[2][-1] == 'F':
        moves.append(f"Move F from 3 to 1")
        current[2].pop()
        current[0].append('F')
    
    while current[1][-1] == 'K':
        moves.append(f"Move K from 2 to 1")
        current[1].pop()
        current[0].append('K')
        
    # Move blocks to their target positions
    for target_stack_idx, target_stack in enumerate(goal):
        for target_block in target_stack:
            # Find where the block currently is
            for curr_stack_idx, curr_stack in enumerate(current):
                if curr_stack and target_block in curr_stack:
                    # If block is not at top, move blocks above it
                    while curr_stack[-1] != target_block:
                        temp_block = curr_stack[-1]
                        # Find temporary position (prefer empty or non-goal position)
                        temp_stack = next((i for i in range(3) 
                                        if i != curr_stack_idx and i != target_stack_idx), 0)
                        moves.append(f"Move {temp_block} from {curr_stack_idx + 1} to {temp_stack + 1}")
                        curr_stack.pop()
                        current[temp_stack].append(temp_block)
                    
                    # Move the target block
                    if curr_stack_idx != target_stack_idx:
                        moves.append(f"Move {target_block} from {curr_stack_idx + 1} to {target_stack_idx + 1}")
                        curr_stack.pop()
                        current[target_stack_idx].append(target_block)
                    break

print('\n'.join(find_solution()))