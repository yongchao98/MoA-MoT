def print_solution():
    moves = []
    
    # Step 1: Move C from stack1 to stack2
    moves.append("Move C from 1 to 2")
    
    # Step 2: Move E from stack1 to stack4
    moves.append("Move E from 1 to 4")
    
    # Step 3: Move F from stack3 to stack2
    moves.append("Move F from 3 to 2")
    
    # Step 4: Move B from stack4 to stack3
    moves.append("Move B from 4 to 3")
    
    # Step 5: Move E from stack4 to stack3
    moves.append("Move E from 4 to 3")
    
    # Step 6: Move F from stack2 to stack3
    moves.append("Move F from 2 to 3")
    
    print("\n".join(moves))

print_solution()