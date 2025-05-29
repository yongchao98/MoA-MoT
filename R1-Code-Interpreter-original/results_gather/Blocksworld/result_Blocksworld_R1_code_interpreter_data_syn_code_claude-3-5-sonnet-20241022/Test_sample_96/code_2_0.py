def print_state(stacks):
    for i, stack in enumerate(stacks, 1):
        print(f"Stack {i}: {' '.join(stack)}")

def make_move(stacks, from_stack, to_stack):
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    return block

def find_solution():
    # Initial state
    stacks = [
        ['E', 'B', 'J', 'I'],
        ['A', 'F', 'K', 'C'],
        ['D', 'G', 'H']
    ]
    
    moves = []
    
    # Moving to goal state: 
    # stack1: A C F J K
    # stack2: D G
    # stack3: B E H I
    
    # First clear space in stack2
    moves.append(('C', 2, 3))  # Move C from stack2 to stack3
    make_move(stacks, 2, 3)
    
    moves.append(('K', 2, 3))  # Move K from stack2 to stack3
    make_move(stacks, 2, 3)
    
    moves.append(('F', 2, 3))  # Move F from stack2 to stack3
    make_move(stacks, 2, 3)
    
    # Move I from stack1 to stack2
    moves.append(('I', 1, 2))
    make_move(stacks, 1, 2)
    
    # Move J from stack1 to stack2
    moves.append(('J', 1, 2))
    make_move(stacks, 1, 2)
    
    # Move B from stack1 to stack3
    moves.append(('B', 1, 3))
    make_move(stacks, 1, 3)
    
    # Move E from stack1 to stack3
    moves.append(('E', 1, 3))
    make_move(stacks, 1, 3)
    
    # Build stack2
    moves.append(('J', 2, 3))  # Temporarily move J
    make_move(stacks, 2, 3)
    
    moves.append(('I', 2, 3))  # Move I to final position
    make_move(stacks, 2, 3)
    
    # Move H to position
    moves.append(('H', 3, 2))
    make_move(stacks, 3, 2)
    
    # Move G to position
    moves.append(('G', 3, 2))
    make_move(stacks, 3, 2)
    
    # Move D to position
    moves.append(('D', 3, 2))
    make_move(stacks, 3, 2)
    
    # Build stack1
    moves.append(('F', 3, 1))
    make_move(stacks, 3, 1)
    
    moves.append(('K', 3, 1))
    make_move(stacks, 3, 1)
    
    moves.append(('J', 3, 1))
    make_move(stacks, 3, 1)
    
    moves.append(('C', 3, 1))
    make_move(stacks, 3, 1)
    
    moves.append(('A', 2, 1))
    make_move(stacks, 2, 1)
    
    # Format and print the moves
    for block, source, dest in moves:
        print(f"Move {block} from {source} to {dest}")

find_solution()