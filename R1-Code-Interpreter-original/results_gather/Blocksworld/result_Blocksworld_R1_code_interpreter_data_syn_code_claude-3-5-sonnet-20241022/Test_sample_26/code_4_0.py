def print_solution():
    moves = []
    # Clear stack 1 first
    moves.append("Move C from 1 to 3")
    moves.append("Move J from 1 to 2")
    moves.append("Move K from 1 to 3")
    moves.append("Move B from 1 to 2")
    
    # Clear access to D
    moves.append("Move A from 3 to 2")
    moves.append("Move I from 3 to 1")
    moves.append("Move G from 3 to 2")
    moves.append("Move E from 3 to 2")
    
    # Build stack 1 (C D F I K)
    moves.append("Move D from 3 to 1")
    moves.append("Move F from 2 to 1")
    moves.append("Move I from 1 to 3")  # Temporarily move I
    moves.append("Move F from 1 to 3")  # Temporarily move F
    moves.append("Move D from 1 to 3")  # Temporarily move D
    moves.append("Move C from 3 to 1")  # Place C at bottom
    moves.append("Move D from 3 to 1")  # Build up stack 1
    moves.append("Move F from 3 to 1")
    moves.append("Move I from 3 to 1")
    moves.append("Move K from 3 to 1")
    
    # Arrange stack 2 (B H J)
    moves.append("Move E from 2 to 3")
    moves.append("Move G from 2 to 3")
    moves.append("Move A from 2 to 3")
    
    print('\n'.join(moves))

print_solution()