def print_state(stacks):
    for i, stack in enumerate(stacks, 1):
        print(f"stack{i}: {' '.join(stack)}")
    print()

def solve_blocksworld():
    moves = []
    
    # Initial state
    stacks = [
        ['D', 'E'],        # stack1
        ['A', 'I', 'G'],   # stack2
        ['C', 'B'],        # stack3
        ['K', 'F', 'J', 'H']  # stack4
    ]
    
    # Move H from stack4 to stack2
    moves.append(('H', 4, 2))
    
    # Move J from stack4 to stack2
    moves.append(('J', 4, 2))
    
    # Move F from stack4 to stack2
    moves.append(('F', 4, 2))
    
    # Move K from stack4 to stack2
    moves.append(('K', 4, 2))
    
    # Move B from stack3 to stack3 (final position)
    moves.append(('B', 3, 3))
    
    # Move C from stack3 to stack1
    moves.append(('C', 3, 1))
    
    # Move G from stack2 to stack1
    moves.append(('G', 2, 1))
    
    # Move E from stack1 to stack2
    moves.append(('E', 1, 2))
    
    # Move D from stack1 to stack1 (temporary)
    moves.append(('D', 1, 4))
    
    # Move A from stack2 to stack1
    moves.append(('A', 2, 1))
    
    # Move D back to stack1
    moves.append(('D', 4, 1))

    result = ""
    for block, source, dest in moves:
        result += f"Move {block} from {source} to {dest}\n"
    print(result.strip())

solve_blocksworld()