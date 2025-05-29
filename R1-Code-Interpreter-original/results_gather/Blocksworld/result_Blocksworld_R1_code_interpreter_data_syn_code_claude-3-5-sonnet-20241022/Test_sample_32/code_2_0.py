def print_final_solution():
    moves = []
    moves.append("Move C from 1 to 2")
    moves.append("Move E from 2 to 3")
    moves.append("Move B from 1 to 3")
    moves.append("Move D from 1 to 3")
    moves.append("Move G from 1 to 3")
    moves.append("Move A from 2 to 1")
    moves.append("Move B from 3 to 1")
    moves.append("Move D from 3 to 1")
    moves.append("Move E from 3 to 1")
    moves.append("Move G from 3 to 1")
    moves.append("Move F from 2 to 3")
    
    print("<<<" + "\n".join(moves) + ">>>")

print_final_solution()