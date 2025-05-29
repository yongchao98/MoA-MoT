def print_valid_solution():
    # Initial state
    stacks = [
        ['F', 'A'],
        ['K', 'D', 'E', 'B', 'I', 'C'],
        ['H', 'J', 'G']
    ]
    
    moves = []
    
    # Move C to final position
    moves.append("Move C from 2 to 3")
    # Move I to final position
    moves.append("Move I from 2 to 3")
    # Move B to stack2 (final position)
    moves.append("Move B from 2 to 1")
    # Move E to final position
    moves.append("Move E from 2 to 3")
    # Move D to stack2 (final position)
    moves.append("Move D from 2 to 2")
    # Move K to final position
    moves.append("Move K from 2 to 3")
    # Clear stack1
    moves.append("Move A from 1 to 3")
    moves.append("Move F from 1 to 3")
    # Move B to stack2
    moves.append("Move B from 1 to 2")
    # Move G to stack1
    moves.append("Move G from 3 to 1")
    
    print("<<<" + "\n".join(moves) + ">>>")

print_valid_solution()