def print_solution():
    moves = []
    
    # Get B D in position first
    moves.append("Move F from 4 to 3")
    moves.append("Move E from 4 to 3")
    moves.append("Move A from 4 to 2")
    moves.append("Move B from 1 to 4")
    moves.append("Move D from 1 to 3")
    moves.append("Move B from 4 to 1")
    moves.append("Move D from 3 to 1")
    
    # Now arrange A C E F in stack2
    moves.append("Move E from 3 to 2")
    moves.append("Move F from 3 to 2")
    moves.append("Move C from 3 to 2")
    
    print("<<<" + "\n".join(moves) + ">>>")

print_solution()