def print_absolutely_final_solution():
    moves = []
    
    # First clear stack2 except D
    moves.append("Move C from 2 to 3")  # C to final stack
    moves.append("Move I from 2 to 3")  # I to final stack
    moves.append("Move B from 2 to 1")  # Temporarily store B
    moves.append("Move E from 2 to 3")  # E to final stack
    moves.append("Move K from 2 to 3")  # K to final stack
    # Move blocks from stack1 to final positions
    moves.append("Move A from 1 to 3")  # A to final stack
    moves.append("Move F from 1 to 3")  # F to final stack
    # Get B D in correct order
    moves.append("Move D from 2 to 1")  # Temporarily move D
    moves.append("Move B from 1 to 2")  # B to final position
    moves.append("Move D from 1 to 2")  # D on top of B
    # Get G to its final position
    moves.append("Move G from 3 to 1")  # G to final position
    # Reorder stack3 (H J already at bottom)
    moves.append("Move K from 3 to 1")  # Temporarily store K
    moves.append("Move F from 3 to 2")  # Temporarily store F
    moves.append("Move A from 3 to 1")  # Temporarily store A
    moves.append("Move E from 3 to 2")  # Temporarily store E
    moves.append("Move I from 3 to 1")  # Temporarily store I
    moves.append("Move C from 3 to 2")  # Temporarily store C
    # Now rebuild stack3 in correct order
    moves.append("Move C from 2 to 3")
    moves.append("Move I from 1 to 3")
    moves.append("Move E from 2 to 3")
    moves.append("Move F from 2 to 3")
    moves.append("Move A from 1 to 3")
    moves.append("Move K from 1 to 3")
    
    print("<<<" + "\n".join(moves) + ">>>")

print_absolutely_final_solution()