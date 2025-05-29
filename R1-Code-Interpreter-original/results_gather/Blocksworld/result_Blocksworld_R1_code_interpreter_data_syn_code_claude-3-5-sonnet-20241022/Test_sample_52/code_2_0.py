def print_solution():
    # Initial state:
    # stack1: I J F
    # stack2: B G E C
    # stack3: H D A
    
    # Goal state:
    # stack1: B C H I J
    # stack2: D
    # stack3: A E F G
    
    moves = [
        "Move C from 2 to 3",  # Move C off to make way
        "Move E from 2 to 3",  # Move E to its final position in stack 3
        "Move G from 2 to 3",  # Move G to its final position in stack 3
        "Move B from 2 to 1",  # Move B to its final position in stack 1
        "Move F from 1 to 3",  # Move F to its final position in stack 3
        "Move J from 1 to 1",  # Temporarily move J
        "Move I from 1 to 1",  # Temporarily move I
        "Move D from 3 to 2",  # Move D to its final position in stack 2
        "Move H from 3 to 1",  # Move H to its final position in stack 1
        "Move I from 1 to 1",  # Move I back to its final position
        "Move J from 1 to 1"   # Move J back to its final position
    ]
    
    # Filter out unnecessary moves (where source and destination are the same)
    final_moves = [move for move in moves if move.split()[3] != move.split()[5]]
    
    print("<<<" + "\n".join(final_moves) + ">>>")

print_solution()