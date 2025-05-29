def get_top(stack):
    return stack[-1] if stack else None

def is_valid_move(stacks, from_stack, to_stack):
    if not stacks[from_stack-1]:  # Source stack is empty
        return False
    return True

def make_move(stacks, move):
    from_stack, to_stack = move
    if not is_valid_move(stacks, from_stack, to_stack):
        return None
    new_stacks = [stack[:] for stack in stacks]
    block = new_stacks[from_stack-1].pop()
    new_stacks[to_stack-1].append(block)
    return new_stacks, block

def solve_blocksworld():
    # Initial and goal states
    initial = [['F', 'I', 'H'], ['G', 'A', 'E', 'C', 'B'], ['D']]
    goal = [['A', 'C'], ['E', 'F', 'G'], ['B', 'D', 'H', 'I']]
    
    moves = []
    current = [stack[:] for stack in initial]
    
    # Strategy: First clear space in stack1 and stack2 for final configuration
    while len(current[0]) > 2:  # Clear stack1 to have only 2 blocks
        block = current[0].pop()
        current[2].append(block)
        moves.append(f"Move {block} from 1 to 3")
    
    while len(current[1]) > 3:  # Clear stack2 to have only 3 blocks
        block = current[1].pop()
        current[2].append(block)
        moves.append(f"Move {block} from 2 to 3")
    
    # Move blocks to their final positions
    # First handle stack1 (A C)
    while current[0] != goal[0]:
        # If we need A and it's accessible
        if 'A' in current[1] and current[1][-1] == 'A':
            current[1].pop()
            current[0].append('A')
            moves.append("Move A from 2 to 1")
        # If we need C and it's accessible
        elif 'C' in current[1] and current[1][-1] == 'C':
            current[1].pop()
            current[0].append('C')
            moves.append("Move C from 2 to 1")
        else:
            # Move top block to stack3 if it's not needed
            if current[1]:
                block = current[1].pop()
                current[2].append(block)
                moves.append(f"Move {block} from 2 to 3")
    
    # Then handle stack2 (E F G)
    target_stack2 = ['E', 'F', 'G']
    for target in target_stack2:
        # Find the block and move it to stack2
        for i in range(3):
            if current[i] and current[i][-1] == target:
                current[i].pop()
                current[1].append(target)
                moves.append(f"Move {target} from {i+1} to 2")
                break
    
    print("<<<" + "\n".join(moves) + ">>>")

solve_blocksworld()