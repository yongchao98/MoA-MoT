def print_solution():
    # Initial state: stack1: H E D F A, stack2: G B, stack3: C
    # Goal state: stack1: G, stack2: A B D E, stack3: C F H
    
    # Direct solution based on analysis:
    moves = [
        "Move H from 1 to 3",  # Move top blocks to clear access
        "Move E from 1 to 2",
        "Move D from 1 to 2",
        "Move F from 1 to 3",
        "Move A from 1 to 2",  # Move A to its final position
        "Move G from 2 to 1",  # Move G to its final position
        "Move B from 2 to 2"   # B is already in correct stack
    ]
    
    print('\n'.join(moves))

print_solution()