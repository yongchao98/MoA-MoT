def solve_blocksworld():
    moves = []
    
    # Initial state
    stacks = [
        ['D', 'E'],        # stack1
        ['A', 'I', 'G'],   # stack2
        ['C', 'B'],        # stack3
        ['K', 'F', 'J', 'H']  # stack4
    ]
    
    # Clear stack4 first (top to bottom)
    moves.append(('H', 4, 1))  # temporary
    moves.append(('J', 4, 2))  # temporary
    moves.append(('F', 4, 1))  # temporary
    moves.append(('K', 4, 2))  # temporary
    
    # Clear stack3 to get B
    moves.append(('B', 3, 3))  # B to final position
    
    # Move C to stack1
    moves.append(('C', 3, 1))  # C to final position
    
    # Clear stack2 and build stack1
    moves.append(('G', 2, 1))  # temporary
    moves.append(('I', 2, 2))  # I to final position
    moves.append(('A', 2, 1))  # A to final position
    
    # Build stack2
    moves.append(('E', 1, 2))  # E to final position
    moves.append(('F', 1, 2))  # F to final position
    moves.append(('H', 1, 2))  # H to final position
    # I is already in position
    moves.append(('J', 2, 2))  # J to final position
    moves.append(('K', 2, 2))  # K to final position
    
    # Complete stack1
    moves.append(('G', 1, 4))  # temporary
    moves.append(('D', 1, 1))  # D to final position
    moves.append(('G', 4, 1))  # G to final position

    result = ""
    for block, source, dest in moves:
        result += f"Move {block} from {source} to {dest}\n"
    print(result.strip())

solve_blocksworld()