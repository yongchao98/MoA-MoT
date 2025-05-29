def print_solution():
    # Initial state
    stacks = [
        ['J'],                          # stack 1
        ['C', 'K', 'F', 'E', 'B', 'I'], # stack 2
        ['G', 'A'],                     # stack 3
        ['H', 'D']                      # stack 4
    ]
    
    moves = []
    
    # Carefully ordered sequence to reach goal state
    # Step 1: Clear space and prepare for F G H I J in stack1
    moves.append("Move I from 2 to 4")  # Clear I temporarily
    moves.append("Move B from 2 to 3")  # Move B to final position
    moves.append("Move E from 2 to 4")  # Clear E
    moves.append("Move F from 2 to 1")  # Start building stack1
    
    # Step 2: Build stack2 (C D E K)
    moves.append("Move K from 2 to 2")  # Keep K in position
    moves.append("Move D from 4 to 2")  # Move D to stack2
    moves.append("Move E from 4 to 2")  # Move E back to stack2
    
    # Step 3: Build stack1 (F G H I J)
    moves.append("Move G from 3 to 1")  # Continue stack1
    moves.append("Move H from 4 to 1")  # Continue stack1
    moves.append("Move I from 4 to 1")  # Continue stack1
    
    # Step 4: Ensure A B in stack3 (already done by previous moves)
    
    print("<<<" + "\n".join(moves) + ">>>")

print_solution()