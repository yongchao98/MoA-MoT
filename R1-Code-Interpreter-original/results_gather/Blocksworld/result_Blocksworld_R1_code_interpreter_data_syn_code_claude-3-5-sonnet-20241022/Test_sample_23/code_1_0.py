def solve_blocksworld():
    # Initial state
    stacks = [
        ['B', 'I'],      # stack 1
        ['G', 'F'],      # stack 2
        ['C', 'D', 'E'], # stack 3
        ['A', 'H']       # stack 4
    ]
    
    moves = []
    
    # Step 1: Clear space for E in stack1
    if 'I' in stacks[0]:
        moves.append("Move I from 1 to 3")
    
    # Step 2: Move E to stack1
    moves.append("Move E from 3 to 1")
    
    # Step 3: Set up stack2 (A D)
    moves.append("Move D from 3 to 2")
    moves.append("Move F from 2 to 3")
    moves.append("Move A from 4 to 2")
    
    # Step 4: Build stack3 (B C F G H I)
    moves.append("Move H from 4 to 3")
    moves.append("Move G from 2 to 3")
    moves.append("Move B from 1 to 3")
    
    print('\n'.join(moves))

solve_blocksworld()