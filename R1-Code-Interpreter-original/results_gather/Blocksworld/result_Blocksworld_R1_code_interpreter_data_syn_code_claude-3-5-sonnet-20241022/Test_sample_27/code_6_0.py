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
    moves.append(('C', 3, 1))
    
    # Build stack1
    moves.append(('G', 2, 4))
    moves.append(('I', 2, 2))
    moves.append(('A', 2, 1))
    
    # Arrange blocks for stack2
    moves.append(('F', 1, 2))
    moves.append(('H', 1, 2))
    moves.append(('E', 1, 2))
    
    # Move D temporarily
    moves.append(('D', 1, 4))
    
    # Complete stack1
    moves.append(('D', 4, 1))
    moves.append(('G', 4, 1))

    result = ""
    for block, source, dest in moves:
        result += f"Move {block} from {source} to {dest}\n"
    print(result.strip())

solve_blocksworld()