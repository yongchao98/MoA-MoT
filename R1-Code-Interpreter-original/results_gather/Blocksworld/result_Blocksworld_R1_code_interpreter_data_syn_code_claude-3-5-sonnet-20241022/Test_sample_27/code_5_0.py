def solve_blocksworld():
    moves = []
    
    # Initial state
    stacks = [
        ['D', 'E'],        # stack1
        ['A', 'I', 'G'],   # stack2
        ['C', 'B'],        # stack3
        ['K', 'F', 'J', 'H']  # stack4
    ]
    
    # Clear stack4
    moves.append(('H', 4, 1))
    moves.append(('J', 4, 2))
    moves.append(('F', 4, 1))
    moves.append(('K', 4, 2))
    
    # Get B to final position
    moves.append(('C', 3, 1))  # Move C first
    moves.append(('B', 3, 3))  # Now B can move to final position
    
    # Build stack1
    moves.append(('G', 2, 4))  # temporary
    moves.append(('I', 2, 2))  # temporary
    moves.append(('A', 2, 1))  # A to final position
    
    # Move C to position
    moves.append(('F', 1, 2))
    moves.append(('H', 1, 2))
    moves.append(('C', 1, 1))  # C to final position
    
    # Build stack2
    moves.append(('E', 1, 2))  # E to final position
    moves.append(('D', 1, 4))  # temporary
    moves.append(('F', 2, 2))  # F to final position
    moves.append(('H', 2, 2))  # H to final position
    moves.append(('I', 2, 2))  # I to final position
    moves.append(('J', 2, 2))  # J to final position
    moves.append(('K', 2, 2))  # K to final position
    
    # Complete stack1
    moves.append(('D', 4, 1))  # D to final position
    moves.append(('G', 4, 1))  # G to final position

    result = ""
    for block, source, dest in moves:
        result += f"Move {block} from {source} to {dest}\n"
    print(result.strip())

solve_blocksworld()