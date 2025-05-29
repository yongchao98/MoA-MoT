def print_final_solution():
    moves = []
    
    # First get G ready to move
    moves.append("Move C from 2 to 1")
    moves.append("Move I from 2 to 1")
    moves.append("Move B from 2 to 3")
    moves.append("Move E from 2 to 1")
    moves.append("Move K from 2 to 3")
    # Now D is exposed in stack2
    moves.append("Move E from 1 to 3")
    moves.append("Move I from 1 to 3")
    moves.append("Move C from 1 to 3")
    moves.append("Move A from 1 to 3")
    moves.append("Move F from 1 to 3")
    # Move B to stack2 on top of D
    moves.append("Move B from 3 to 2")
    # Finally move G to stack1
    moves.append("Move G from 3 to 1")
    
    print("<<<" + "\n".join(moves) + ">>>")

print_final_solution()