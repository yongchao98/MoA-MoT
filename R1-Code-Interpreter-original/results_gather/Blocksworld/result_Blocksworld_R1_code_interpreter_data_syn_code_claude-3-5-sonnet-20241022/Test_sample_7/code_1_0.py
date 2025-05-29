def print_solution():
    # Initial state: C | GABF | DH | E
    # Goal state: AE | CD | BFGH
    
    moves = []
    
    # Step 1: Move F, B to stack 3 to prepare for G
    moves.append("Move F from 2 to 3")
    moves.append("Move B from 2 to 3")
    
    # Step 2: Move A to stack 1 for final position
    moves.append("Move A from 2 to 1")
    
    # Step 3: Move G to stack 3 for final position
    moves.append("Move G from 2 to 3")
    
    # Step 4: Move C to stack 2 for final position
    moves.append("Move C from 1 to 2")
    
    # Step 5: Move D to stack 2 for final position
    moves.append("Move D from 3 to 2")
    
    # Step 6: Move E to stack 1 for final position
    moves.append("Move E from 4 to 1")
    
    print("<<<" + "\n".join(moves) + ">>>")

print_solution()