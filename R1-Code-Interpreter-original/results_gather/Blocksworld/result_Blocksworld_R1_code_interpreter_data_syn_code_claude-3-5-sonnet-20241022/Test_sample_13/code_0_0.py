def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack is not empty
    if not stacks[from_stack-1]:
        return False
    return True

def make_move(stacks, from_stack, to_stack):
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    return block

def print_state(stacks):
    for i, stack in enumerate(stacks, 1):
        print(f"Stack {i}: {' '.join(stack)}")
    print()

def find_solution():
    # Initial state
    stacks = [
        ['K', 'I', 'G'],
        ['A', 'F', 'B'],
        ['J', 'E'],
        ['C', 'H', 'D']
    ]
    
    # Goal state
    goal = [
        ['B', 'E', 'F', 'H'],
        ['A', 'D', 'I'],
        ['C', 'G', 'J', 'K']
    ]
    
    moves = []
    
    # Helper function to get top block
    def get_top(stack_idx):
        return stacks[stack_idx-1][-1] if stacks[stack_idx-1] else None
    
    # Strategy: First clear space by moving blocks to appropriate positions
    while stacks != goal:
        # Move D from stack4 to stack2 (it belongs there)
        if 'D' in stacks[3] and get_top(4) == 'D':
            moves.append(f"Move D from 4 to 2")
            make_move(stacks, 4, 2)
            continue
            
        # Move H from stack4 to stack1 (it belongs there)
        if 'H' in stacks[3] and get_top(4) == 'H':
            moves.append(f"Move H from 4 to 1")
            make_move(stacks, 4, 1)
            continue
            
        # Move G from stack1 to stack3
        if 'G' in stacks[0] and get_top(1) == 'G':
            moves.append(f"Move G from 1 to 3")
            make_move(stacks, 1, 3)
            continue
            
        # Move I from stack1 to stack2
        if 'I' in stacks[0] and get_top(1) == 'I':
            moves.append(f"Move I from 1 to 2")
            make_move(stacks, 1, 2)
            continue
            
        # Move K from stack1 to stack3
        if 'K' in stacks[0] and get_top(1) == 'K':
            moves.append(f"Move K from 1 to 3")
            make_move(stacks, 1, 3)
            continue
            
        # Move B from stack2 to stack1
        if 'B' in stacks[1] and get_top(2) == 'B':
            moves.append(f"Move B from 2 to 1")
            make_move(stacks, 2, 1)
            continue
            
        # Move F from stack2 to stack1
        if 'F' in stacks[1] and get_top(2) == 'F':
            moves.append(f"Move F from 2 to 1")
            make_move(stacks, 2, 1)
            continue
            
        # Move E from stack3 to stack1
        if 'E' in stacks[2] and get_top(3) == 'E':
            moves.append(f"Move E from 3 to 1")
            make_move(stacks, 3, 1)
            continue
            
        # Move J from stack3 to stack3
        if 'J' in stacks[2] and get_top(3) == 'J':
            moves.append(f"Move J from 3 to 3")
            make_move(stacks, 3, 3)
            continue
            
        # Move C from stack4 to stack3
        if 'C' in stacks[3] and get_top(4) == 'C':
            moves.append(f"Move C from 4 to 3")
            make_move(stacks, 4, 3)
            continue
            
    print("\n".join(moves))

find_solution()