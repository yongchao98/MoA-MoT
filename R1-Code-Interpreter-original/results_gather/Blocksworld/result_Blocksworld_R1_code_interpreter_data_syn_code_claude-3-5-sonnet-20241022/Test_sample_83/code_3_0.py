def print_moves():
    # Initial state:
    # stack1: G A E C
    # stack2: H I F
    # stack3: K
    # stack4: J B D
    
    # Goal state:
    # stack1: A C
    # stack2: J
    # stack3: B D E F G H I K
    
    moves = [
        # First clear stack 4
        "Move D from 4 to 3",
        "Move B from 4 to 3",
        "Move J from 4 to 2",
        
        # Move blocks from stack 1
        "Move C from 1 to 4",  # temporarily store C
        "Move E from 1 to 3",
        "Move A from 1 to 4",  # temporarily store A
        "Move G from 1 to 3",
        
        # Move blocks from stack 2
        "Move F from 2 to 3",
        "Move I from 2 to 3",
        "Move H from 2 to 3",
        
        # Move K from stack 3 to complete the sequence
        "Move K from 3 to 3",  # reposition K
        
        # Finally, build stack 1
        "Move A from 4 to 1",
        "Move C from 4 to 1"
    ]
    
    print("<<<")
    print("\n".join(moves))
    print(">>>")

print_moves()