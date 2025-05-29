def print_solution():
    moves = []
    
    # Clear space and get B D in position
    moves.append("Move F from 4 to 3")
    moves.append("Move E from 4 to 2")
    moves.append("Move A from 4 to 2")
    moves.append("Move B from 1 to 4")
    moves.append("Move D from 1 to 3")
    
    # Get A C E F in position
    moves.append("Move A from 2 to 1")
    moves.append("Move E from 2 to 3")
    moves.append("Move A from 1 to 2")
    moves.append("Move E from 3 to 2")
    moves.append("Move F from 3 to 2")
    moves.append("Move B from 4 to 1")
    moves.append("Move D from 3 to 1")
    moves.append("Move C from 3 to 2")
    
    print("<<<" + "\n".join(moves) + ">>>")

print_solution()