def print_solution():
    moves = []
    
    # Build stack1 (F I)
    moves.append("Move I from 3 to 4")
    moves.append("Move E from 3 to 4")
    moves.append("Move C from 3 to 4")
    moves.append("Move F from 3 to 1")
    moves.append("Move I from 4 to 1")
    
    # Build stack2 (A B C D G)
    moves.append("Move J from 2 to 4")  # Move J temporarily
    moves.append("Move D from 1 to 4")
    moves.append("Move B from 4 to 2")
    moves.append("Move C from 4 to 2")
    moves.append("Move D from 4 to 2")
    moves.append("Move G from 4 to 2")
    
    # Build stack3 (E H J) - in correct order
    moves.append("Move E from 4 to 3")  # Place E first
    moves.append("Move J from 4 to 3")  # Place J last
    
    print("<<<")
    print("\n".join(moves))
    print(">>>")

print_solution()