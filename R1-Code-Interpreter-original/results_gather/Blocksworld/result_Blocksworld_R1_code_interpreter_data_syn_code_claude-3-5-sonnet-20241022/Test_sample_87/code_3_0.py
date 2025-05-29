def print_solution():
    moves = []
    
    # Step 1: Move F from stack4 (top)
    moves.append("Move F from 4 to 3")
    
    # Step 2: Move E from stack4 (top)
    moves.append("Move E from 4 to 2")
    
    # Step 3: Move A from stack4 (top)
    moves.append("Move A from 4 to 1")
    
    # Step 4: Move B from stack1 (top)
    moves.append("Move B from 1 to 2")
    
    # Step 5: Move D from stack1 (top)
    moves.append("Move D from 1 to 4")
    
    # Step 6: Move A from stack1 (top)
    moves.append("Move A from 1 to 2")
    
    # Step 7: Move F from stack3 (top)
    moves.append("Move F from 3 to 2")
    
    # Step 8: Move B from stack2 (top)
    moves.append("Move B from 2 to 1")
    
    # Step 9: Move D from stack4 (top)
    moves.append("Move D from 4 to 1")
    
    # Step 10: Move C from stack3 (top)
    moves.append("Move C from 3 to 2")
    
    print("<<<" + "\n".join(moves) + ">>>")

print_solution()