def solve_blocksworld():
    moves = []
    
    # Initial state
    stacks = [
        ['D', 'E'],        # stack1
        ['A', 'I', 'G'],   # stack2
        ['C', 'B'],        # stack3
        ['K', 'F', 'J', 'H']  # stack4
    ]
    
    # First clear stack4 in correct order
    moves.append(('H', 4, 2))  # temporary storage
    moves.append(('J', 4, 1))  # temporary storage
    moves.append(('F', 4, 2))  # will be in final position
    moves.append(('K', 4, 2))  # will be in final position
    
    # Move B to its final position
    moves.append(('B', 3, 3))
    
    # Build stack1
    moves.append(('C', 3, 1))  # C to final position
    moves.append(('G', 2, 1))  # G to final position
    moves.append(('E', 1, 2))  # Move E to stack2 (final position)
    moves.append(('D', 1, 1))  # D to final position
    moves.append(('J', 1, 2))  # J to stack2 (final position)
    moves.append(('I', 2, 2))  # I to correct position
    moves.append(('H', 2, 2))  # H to correct position
    moves.append(('A', 2, 1))  # A to final position

    result = ""
    for block, source, dest in moves:
        result += f"Move {block} from {source} to {dest}\n"
    print(result.strip())

solve_blocksworld()