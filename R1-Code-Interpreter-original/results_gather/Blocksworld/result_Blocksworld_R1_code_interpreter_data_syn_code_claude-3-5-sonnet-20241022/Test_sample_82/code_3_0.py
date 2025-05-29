def print_solution():
    moves = []
    
    # First clear space in stack2 for B D
    moves.append("Move C from 2 to 3")
    moves.append("Move I from 2 to 3")
    moves.append("Move B from 2 to 1")
    moves.append("Move E from 2 to 3")
    # K and D are in stack2, need to move K
    moves.append("Move K from 2 to 3")
    # Now clear stack1 except B
    moves.append("Move A from 1 to 3")
    moves.append("Move F from 1 to 3")
    # Move B to its final position
    moves.append("Move B from 1 to 2")
    # Finally move G to stack1
    moves.append("Move G from 3 to 1")
    
    print("<<<" + "\n".join(moves) + ">>>")

print_solution()