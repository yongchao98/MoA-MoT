def print_solution():
    # Initial state
    stacks = [
        ['B', 'K', 'J', 'C'],  # stack 1
        ['F', 'H'],           # stack 2
        ['D', 'E', 'G', 'I', 'A']  # stack 3
    ]
    
    # Goal state
    goal = [
        ['C', 'D', 'F', 'I', 'K'],  # stack 1
        ['B', 'H', 'J'],           # stack 2
        ['A', 'E', 'G']            # stack 3
    ]
    
    moves = []
    
    # Direct solution based on analysis
    # First clear space for C, D, F, I, K in stack 1
    moves.append("Move C from 1 to 3")  # Save C temporarily
    moves.append("Move J from 1 to 2")  # Move J to its goal stack
    moves.append("Move K from 1 to 3")  # Save K temporarily
    moves.append("Move B from 1 to 2")  # Move B to its goal stack
    
    # Now build stack 1
    moves.append("Move A from 3 to 1")  # Clear access to I
    moves.append("Move I from 3 to 1")  # Move I to stack 1
    moves.append("Move G from 3 to 3")  # Clear access to E
    moves.append("Move E from 3 to 3")  # Clear access to D
    moves.append("Move D from 3 to 1")  # Move D to stack 1
    moves.append("Move K from 3 to 1")  # Move K to stack 1
    
    # Complete stack 2
    moves.append("Move H from 2 to 3")  # Temporarily move H
    moves.append("Move J from 2 to 2")  # Keep J in place
    moves.append("Move B from 2 to 2")  # Keep B in place
    moves.append("Move F from 2 to 1")  # Move F to stack 1
    
    print('\n'.join(moves))

print_solution()