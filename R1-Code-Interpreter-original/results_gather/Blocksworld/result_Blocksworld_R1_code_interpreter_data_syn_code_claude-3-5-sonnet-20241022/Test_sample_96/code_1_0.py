def print_state(stacks):
    for i, stack in enumerate(stacks, 1):
        print(f"Stack {i}: {' '.join(stack)}")

def get_top_block(stack):
    return stack[-1] if stack else None

def is_valid_move(stacks, from_stack, to_stack):
    if not stacks[from_stack-1]:  # Check if source stack is empty
        return False
    return True

def make_move(stacks, from_stack, to_stack):
    if not is_valid_move(stacks, from_stack, to_stack):
        return None
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
    
    # Step by step solution
    # First move I from stack1 to stack3
    moves.append(('I', 1, 3))
    make_move(stacks, 1, 3)
    
    # Move J from stack1 to stack2
    moves.append(('J', 1, 2))
    make_move(stacks, 1, 2)
    
    # Move C from stack2 to stack3
    moves.append(('C', 2, 3))
    make_move(stacks, 2, 3)
    
    # Move K from stack2 to stack1
    moves.append(('K', 2, 1))
    make_move(stacks, 2, 1)
    
    # Move F from stack2 to stack1
    moves.append(('F', 2, 1))
    make_move(stacks, 2, 1)
    
    # Move A from stack2 to stack1
    moves.append(('A', 2, 1))
    make_move(stacks, 2, 1)
    
    # Move H from stack3 to stack2
    moves.append(('H', 3, 2))
    make_move(stacks, 3, 2)
    
    # Move C from stack3 to stack1
    moves.append(('C', 3, 1))
    make_move(stacks, 3, 1)
    
    # Move G from stack3 to stack2
    moves.append(('G', 3, 2))
    make_move(stacks, 3, 2)
    
    # Move D from stack3 to stack2
    moves.append(('D', 3, 2))
    make_move(stacks, 3, 2)
    
    # Move B from stack1 to stack3
    moves.append(('B', 1, 3))
    make_move(stacks, 1, 3)
    
    # Move E from stack1 to stack3
    moves.append(('E', 1, 3))
    make_move(stacks, 1, 3)
    
    # Format and print the moves
    for block, source, dest in moves:
        print(f"Move {block} from {source} to {dest}")

find_solution()