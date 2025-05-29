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
        # Step 1: Clear stack 4
        (4, 2),  # Move E from 4 to 2 (will be in final position)
        (4, 1),  # Move A from 4 to 1 (temporary)
        (4, 2),  # Move J from 4 to 2 (will be in final position)
        
        # Step 2: Move blocks from stack 1
        (1, 3),  # Move F from 1 to 3 (temporary)
        (1, 3),  # Move B from 1 to 3 (will be in final position)
        (1, 3),  # Move H from 1 to 3 (will be in final position)
        
        # Step 3: Move blocks from stack 3 to stack 1
        (3, 1),  # Move D from 3 to 1
        (3, 1),  # Move C from 3 to 1
        (3, 1),  # Move G from 3 to 1
        
        # Step 4: Move I to final position
        (2, 3),  # Move I from 2 to 3
        
        # Step 5: Move A to final position in stack 1
        (1, 2),  # Move A from 1 to 2 (temporary)
        (1, 3),  # Clear stack 1
        (2, 1),  # Move A to stack 1 (final position)
        
        # Step 6: Rebuild stack 1
        (3, 1),  # Move G back to 1
        (3, 1),  # Move C back to 1
        (3, 1),  # Move D back to 1
        (3, 1),  # Move F back to 1
    ]
    
    solution = []
    current_state = []
    for from_stack, to_stack in moves:
        if is_valid_move(stacks, from_stack, to_stack):
            block = make_move(stacks, from_stack, to_stack)
            if block:
                solution.append(f"Move {block} from {from_stack} to {to_stack}")
                current_state = [stack.copy() for stack in stacks]
    
    print("\n".join(solution))

find_solution()