def print_state(stacks):
    for i, stack in enumerate(stacks, 1):
        print(f"Stack {i}: {' '.join(stack)}")
    print()

def is_valid_move(stacks, from_stack, to_stack):
    if from_stack < 1 or from_stack > len(stacks) or to_stack < 1 or to_stack > len(stacks):
        return False
    return len(stacks[from_stack-1]) > 0

def make_move(stacks, from_stack, to_stack):
    if not is_valid_move(stacks, from_stack, to_stack):
        return None
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    return block

def find_solution():
    # Initial state
    stacks = [
        ['H', 'B', 'F'],    # stack 1
        ['I'],              # stack 2
        ['G', 'C', 'D'],    # stack 3
        ['J', 'A', 'E']     # stack 4
    ]
    
    moves = [
        # First get E and J to stack 2
        (4, 2),  # Move E from 4 to 2
        (4, 1),  # Move A from 4 to 1
        (4, 2),  # Move J from 4 to 2
        
        # Clear stack 1 temporarily
        (1, 4),  # Move F from 1 to 4
        (1, 4),  # Move B from 1 to 4
        (1, 4),  # Move H from 1 to 4
        
        # Build stack 1 from bottom up
        (3, 1),  # Move D from 3 to 1
        (3, 1),  # Move C from 3 to 1
        (3, 1),  # Move G from 3 to 1
        (4, 1),  # Move F from 4 to 1
        (2, 1),  # Move A from 2 to 1
        
        # Build stack 3
        (4, 3),  # Move H from 4 to 3
        (4, 3),  # Move B from 4 to 3
        (2, 3),  # Move I from 2 to 3
    ]
    
    solution = []
    for from_stack, to_stack in moves:
        if is_valid_move(stacks, from_stack, to_stack):
            block = make_move(stacks, from_stack, to_stack)
            if block:
                solution.append(f"Move {block} from {from_stack} to {to_stack}")
    
    print("\n".join(solution))

find_solution()