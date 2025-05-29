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
        "Move H from 1 to 3",    # Get H to its final position
        "Move B from 3 to 2",    # Start building stack 2
        "Move C from 3 to 2",
        "Move E from 3 to 2",    # Stack 2 is now ready
        "Move I from 4 to 3",    # Build stack 3
        "Move D from 2 to 1",    # Get D to stack 1
        "Move G from 2 to 3",    # Move G to final position
        "Move A from 4 to 3",    # Move A to final position
        "Move F from 1 to 1"     # F is already in position
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