def solve_blocks():
    moves = []
    
    # We know the exact sequence needed to reach the goal state:
    # From initial state:
    # stack1: 
    # stack2: B A C
    # stack3: E D
    
    # To goal state:
    # stack1: B E
    # stack2: C
    # stack3: A D
    
    moves = [
        "Move C from 2 to 3",  # stack2: B A, stack3: E D C
        "Move A from 2 to 3",  # stack2: B, stack3: E D C A
        "Move B from 2 to 1",  # stack1: B, stack2: empty, stack3: E D C A
        "Move A from 3 to 3",  # stack1: B, stack2: empty, stack3: E D C A
        "Move C from 3 to 2",  # stack1: B, stack2: C, stack3: E D A
        "Move A from 3 to 3",  # stack1: B, stack2: C, stack3: E D A
        "Move E from 3 to 1"   # stack1: B E, stack2: C, stack3: D A
    ]
    
    print("\n".join(moves))

solve_blocks()