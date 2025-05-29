def print_state(stacks):
    for i, stack in enumerate(stacks):
        print(f"Stack {i+1}: {stack}")
    print()

def is_valid_move(stacks, from_stack, to_stack):
    if from_stack < 1 or from_stack > len(stacks) or to_stack < 1 or to_stack > len(stacks):
        return False
    if not stacks[from_stack-1]:  # Empty source stack
        return False
    return True

def make_move(stacks, from_stack, to_stack):
    if not is_valid_move(stacks, from_stack, to_stack):
        return None
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    return f"Move {block} from {from_stack} to {to_stack}"

def solve_blocksworld():
    # Initial state
    current_state = [
        ['A', 'G'],        # Stack 1
        ['F', 'E'],        # Stack 2
        ['C'],            # Stack 3
        ['H', 'J', 'I', 'B', 'D']  # Stack 4
    ]
    
    # Goal state
    goal_state = [
        ['A', 'B'],        # Stack 1
        ['C', 'E', 'G', 'J'],  # Stack 2
        ['D', 'F', 'H', 'I']   # Stack 3
    ]
    
    moves = []
    
    # Manually planned sequence
    move_sequence = [
        (4, 3),  # Move D from 4 to 3
        (4, 1),  # Move B from 4 to 1
        (4, 3),  # Move I from 4 to 3
        (4, 2),  # Move J from 4 to 2
        (4, 3),  # Move H from 4 to 3
        (1, 2),  # Move G from 1 to 2
        (2, 3),  # Move F from 2 to 3
        (3, 2),  # Move C from 3 to 2
    ]
    
    # Execute moves
    for from_stack, to_stack in move_sequence:
        move = make_move(current_state, from_stack, to_stack)
        if move:
            moves.append(move)
    
    print('\n'.join(moves))

solve_blocksworld()