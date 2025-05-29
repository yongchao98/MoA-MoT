def print_solution():
    moves = []
    
    # Build stack1 (F I)
    moves.append("Move I from 3 to 4")
    moves.append("Move E from 3 to 4")
    moves.append("Move C from 3 to 4")
    moves.append("Move F from 3 to 1")
    moves.append("Move I from 4 to 1")
    
    # Build stack2 (A B C D G) - in correct order
    moves.append("Move J from 2 to 3")  # Clear J from stack2
    moves.append("Move D from 1 to 4")  # Temporarily move D
    moves.append("Move B from 4 to 2")  # Place B on A
    moves.append("Move C from 4 to 2")  # Place C
    moves.append("Move D from 4 to 2")  # Place D
    moves.append("Move G from 4 to 2")  # Place G
    
    # Build stack3 (E H J)
    moves.append("Move E from 4 to 3")
    
    print("<<<")
    print("\n".join(moves))
    print(">>>")

print_solution()