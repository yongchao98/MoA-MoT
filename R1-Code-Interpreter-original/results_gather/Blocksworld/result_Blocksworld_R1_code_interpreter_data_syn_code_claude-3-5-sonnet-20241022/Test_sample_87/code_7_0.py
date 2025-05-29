def print_solution():
    moves = []
    
    # First get B D in position
    moves.append("Move F from 4 to 2")
    moves.append("Move E from 4 to 2")
    moves.append("Move A from 4 to 3")  # temporarily store A
    moves.append("Move B from 1 to 2")  # temporarily store B
    moves.append("Move D from 1 to 4")  # temporarily store D
    
    # Now build stack2 in correct order (A C E F)
    moves.append("Move B from 2 to 3")
    moves.append("Move E from 2 to 4")
    moves.append("Move F from 2 to 4")
    moves.append("Move A from 3 to 2")  # A first
    moves.append("Move C from 3 to 2")  # C second
    moves.append("Move E from 4 to 2")  # E third
    moves.append("Move F from 4 to 2")  # F fourth
    
    # Finally position B D in stack1
    moves.append("Move B from 3 to 1")
    moves.append("Move D from 4 to 1")
    
    print("<<<" + "\n".join(moves) + ">>>")

print_solution()