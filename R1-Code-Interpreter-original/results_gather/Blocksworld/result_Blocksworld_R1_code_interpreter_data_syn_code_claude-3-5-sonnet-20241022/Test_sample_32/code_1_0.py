def print_solution():
    moves = []
    # First clear space for A to be bottom
    moves.append("Move C from 1 to 2")  # stack1: G D B, stack2: F A E C, stack3: 
    moves.append("Move E from 2 to 3")  # stack1: G D B, stack2: F A C, stack3: E
    moves.append("Move B from 1 to 3")  # stack1: G D, stack2: F A C, stack3: E B
    moves.append("Move D from 1 to 3")  # stack1: G, stack2: F A C, stack3: E B D
    moves.append("Move G from 1 to 3")  # stack1: empty, stack2: F A C, stack3: E B D G
    moves.append("Move A from 2 to 1")  # stack1: A, stack2: F C, stack3: E B D G
    moves.append("Move B from 3 to 1")  # stack1: A B, stack2: F C, stack3: E D G
    moves.append("Move D from 3 to 1")  # stack1: A B D, stack2: F C, stack3: E G
    moves.append("Move E from 3 to 1")  # stack1: A B D E, stack2: F C, stack3: G
    moves.append("Move G from 3 to 1")  # stack1: A B D E G, stack2: F C, stack3: empty
    moves.append("Move C from 2 to 2")  # stack1: A B D E G, stack2: C, stack3: F
    moves.append("Move F from 2 to 3")  # Final state achieved
    
    print("<<<" + "\n".join(moves) + ">>>")

print_solution()