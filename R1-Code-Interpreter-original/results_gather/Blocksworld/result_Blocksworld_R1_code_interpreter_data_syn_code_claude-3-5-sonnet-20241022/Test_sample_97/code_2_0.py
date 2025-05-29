def print_moves():
    # Initial state:
    # stack1: H J A B
    # stack2: F I K C
    # stack3: 
    # stack4: G E D
    
    # Goal state:
    # stack1: B C D E I K
    # stack2: A H J
    # stack3: F G
    
    moves = [
        "Move B from 1 to 3",    # Remove B to make A accessible
        "Move A from 1 to 2",    # Move A to its final position in stack2
        "Move J from 1 to 4",    # Temporarily store J
        "Move H from 1 to 2",    # Move H to its final position in stack2
        "Move J from 4 to 2",    # Move J to its final position in stack2
        "Move C from 2 to 1",    # Start building stack1
        "Move K from 2 to 4",    # Temporarily move K
        "Move I from 2 to 1",    # Move I to stack1
        "Move F from 2 to 3",    # Move F to its final position in stack3
        "Move K from 4 to 1",    # Move K to its final position in stack1
        "Move D from 4 to 1",    # Move D to stack1
        "Move E from 4 to 1",    # Move E to stack1
        "Move G from 4 to 3",    # Move G to its final position in stack3
        "Move B from 3 to 1"     # Finally, move B to complete stack1
    ]
    
    for move in moves:
        print(move)

print_moves()