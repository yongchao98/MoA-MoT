def print_solution():
    moves = []
    
    # First clear stack4 and get blocks in position
    moves.append("Move F from 4 to 2")  # F is on top of stack4
    moves.append("Move E from 4 to 3")  # E is now on top of stack4
    moves.append("Move A from 4 to 2")  # A is now on top of stack4
    
    # Now arrange stack1
    moves.append("Move B from 1 to 3")  # B is on top of stack1
    moves.append("Move D from 1 to 4")  # D is now on top of stack1
    
    # Build stack2 in correct order (A C E F)
    moves.append("Move F from 2 to 4")  # F was on top of stack2
    moves.append("Move A from 2 to 3")  # A was on top of stack2
    moves.append("Move E from 3 to 2")  # E was on top of stack3
    moves.append("Move C from 3 to 2")  # C was on top of stack3
    moves.append("Move A from 3 to 2")  # A was on top of stack3
    
    # Finally position B D in stack1
    moves.append("Move B from 3 to 1")  # B was on top of stack3
    moves.append("Move D from 4 to 1")  # D was on top of stack4
    moves.append("Move F from 4 to 2")  # F was on top of stack4
    
    print("<<<" + "\n".join(moves) + ">>>")

print_solution()