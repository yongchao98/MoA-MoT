def print_final_correct_solution():
    moves = []
    
    # Clear stack2 except D
    moves.append("Move C from 2 to 3")
    moves.append("Move I from 2 to 3")
    moves.append("Move B from 2 to 1")
    moves.append("Move E from 2 to 3")
    moves.append("Move K from 2 to 3")
    # Clear stack1 for G
    moves.append("Move A from 1 to 3")
    moves.append("Move F from 1 to 3")
    # Get B and D in correct order
    moves.append("Move D from 2 to 3")
    moves.append("Move B from 1 to 2")
    moves.append("Move D from 3 to 2")
    # Finally place G
    moves.append("Move G from 3 to 1")
    
    print("<<<" + "\n".join(moves) + ">>>")

print_final_correct_solution()