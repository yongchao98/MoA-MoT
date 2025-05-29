def print_state(stacks):
    for i, stack in enumerate(stacks, 1):
        print(f"Stack {i}: {' '.join(stack)}")
    print()

def verify_moves():
    # Initial state
    stacks = [
        ['F', 'H'],      # Stack 1
        ['G', 'D'],      # Stack 2
        ['E', 'C', 'B'], # Stack 3
        ['A', 'I']       # Stack 4
    ]
    
    moves = [
        "Move H from 1 to 4",    # Clear stack 1
        "Move D from 2 to 1",    # Get D to its final position
        "Move H from 4 to 3",    # Move H to final position
        "Move I from 4 to 3",    # Move I to final position
        "Move B from 3 to 2",    # Start building stack 2
        "Move C from 3 to 2",    # Continue building stack 2
        "Move E from 3 to 2",    # Complete stack 2
        "Move G from 2 to 3",    # Move G to stack 3
        "Move A from 4 to 3",    # Complete stack 3
        "Move F from 1 to 2",    # Move F temporarily
        "Move D from 1 to 1",    # Adjust D
        "Move F from 2 to 1"     # Complete stack 1
    ]
    
    print("Initial state:")
    print_state(stacks)
    
    for move in moves:
        # Parse move
        _, block, _, source, _, dest = move.split()
        source = int(source)
        dest = int(dest)
        
        # Verify move is valid
        if not stacks[source-1] or stacks[source-1][-1] != block:
            print(f"Invalid move: {move}")
            print(f"Current state:")
            print_state(stacks)
            return False
            
        # Make move
        stacks[source-1].pop()
        stacks[dest-1].append(block)
        
        print(f"After {move}:")
        print_state(stacks)

verify_moves()