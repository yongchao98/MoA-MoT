def print_solution():
    # Initial state
    moves = []
    
    # Step 1: Move C from stack1 to stack2
    moves.append("Move C from 1 to 2")
    
    # Step 2: Move E from stack1 to stack3
    moves.append("Move E from 1 to 3")
    
    # Step 3: Move B from stack4 to stack2 (temporary)
    moves.append("Move B from 4 to 2")
    
    # Step 4: Move E from stack3 to stack4
    moves.append("Move E from 3 to 4")
    
    # Step 5: Move B from stack2 to stack3
    moves.append("Move B from 2 to 3")
    
    # Step 6: Move E from stack4 to stack3
    moves.append("Move E from 4 to 3")
    
    print("\n".join(moves))

print_solution()