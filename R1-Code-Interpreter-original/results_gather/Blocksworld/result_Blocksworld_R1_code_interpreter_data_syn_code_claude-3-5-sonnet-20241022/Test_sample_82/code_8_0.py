def print_correct_solution():
    moves = []
    
    # First get B and D in position in stack2
    moves.append("Move C from 2 to 3")
    moves.append("Move I from 2 to 3")
    moves.append("Move B from 2 to 1")  # Temporarily store B
    moves.append("Move E from 2 to 3")
    moves.append("Move K from 2 to 3")
    # Now only D remains in stack2
    
    # Clear stack1 for G
    moves.append("Move A from 1 to 3")
    moves.append("Move F from 1 to 3")
    
    # Put B and D in correct order
    moves.append("Move D from 2 to 3")  # Move D temporarily
    moves.append("Move B from 1 to 2")  # Put B in position
    moves.append("Move D from 3 to 2")  # Put D on top of B
    
    # Finally get G to stack1
    moves.append("Move G from 3 to 1")
    
    print("<<<" + "\n".join(moves) + ">>>")

print_correct_solution()