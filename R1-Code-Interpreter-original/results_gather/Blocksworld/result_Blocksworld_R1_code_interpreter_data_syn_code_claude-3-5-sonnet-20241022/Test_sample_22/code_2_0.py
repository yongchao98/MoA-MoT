def print_moves():
    # Initial state:
    # stack1: I A G
    # stack2: D F E
    # stack3: H B C
    
    # Goal state:
    # stack1: B C H
    # stack2: F G I
    # stack3: A D E
    
    moves = [
        "Move G from 1 to 2",  # Remove G from stack1
        "Move C from 3 to 1",  # Start building stack1 goal
        "Move B from 3 to 1",  # Continue building stack1 goal
        "Move E from 2 to 3",  # Move E to its final position
        "Move H from 3 to 2",  # Move H to temporary position
        "Move I from 1 to 2",  # Move I to its final position
        "Move A from 1 to 3",  # Move A to its final position
        "Move D from 2 to 3"   # Move D to its final position
    ]
    
    print("<<<")
    for move in moves:
        print(move)
    print(">>>")

print_moves()