def print_solution():
    # Initial state:
    # stack1: G E B
    # stack2: A
    # stack3: C D H F
    
    # Goal state:
    # stack1: D F G
    # stack2: A C E H
    # stack3: B
    
    # Direct solution with minimal moves
    moves = [
        "Move B from 1 to 3",    # Free up stack1 for D and F
        "Move E from 1 to 2",    # Move E to its goal position in stack2
        "Move F from 3 to 1",    # Move F to its goal position in stack1
        "Move H from 3 to 2",    # Move H to its goal position in stack2
        "Move D from 3 to 1",    # Move D to its goal position in stack1
        "Move G from 1 to 1",    # G is already in stack1, but needs to be on top of F
        "Move C from 3 to 2"     # Move C to its goal position in stack2
    ]
    
    print("<<<" + "\n".join(moves) + ">>>")

print_solution()