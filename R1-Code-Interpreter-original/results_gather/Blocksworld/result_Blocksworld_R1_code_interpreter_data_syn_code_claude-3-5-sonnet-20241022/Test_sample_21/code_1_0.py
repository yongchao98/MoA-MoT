def print_solution():
    # Initial state
    stacks = [
        [],                     # stack1
        ['F'],                 # stack2
        ['G', 'C'],           # stack3
        ['A', 'B', 'D', 'E', 'H']  # stack4
    ]
    
    moves = []
    
    # Move H from top of stack4 to stack2 (will be part of final G,H)
    moves.append("Move H from 4 to 2")
    
    # Move E from stack4 to stack1 (will be top of B,D,E)
    moves.append("Move E from 4 to 1")
    
    # Move D from stack4 to stack1 (will be middle of B,D,E)
    moves.append("Move D from 4 to 1")
    
    # Move B from stack4 to stack1 (will be bottom of B,D,E)
    moves.append("Move B from 4 to 1")
    
    # Move A from stack4 to stack3 (will be bottom of A,C,F)
    moves.append("Move A from 4 to 3")
    
    # Move F from stack2 to stack3 (will be top of A,C,F)
    moves.append("Move F from 2 to 3")
    
    # Print the solution
    print('\n'.join(moves))

print_solution()