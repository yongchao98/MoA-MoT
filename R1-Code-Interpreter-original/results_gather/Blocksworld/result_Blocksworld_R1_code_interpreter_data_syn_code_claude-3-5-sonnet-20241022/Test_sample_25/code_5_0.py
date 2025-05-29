def print_state(stacks):
    for i, stack in enumerate(stacks):
        print(f"Stack {i+1}: {' '.join(stack)}")

def get_top_block(stack):
    return stack[-1] if stack else None

def execute_moves():
    # Initialize stacks as lists to maintain state
    stacks = [
        ['A', 'G'],        # Stack 1
        ['F', 'E'],        # Stack 2
        ['C'],             # Stack 3
        ['H', 'J', 'I', 'B', 'D']  # Stack 4
    ]
    
    # Hardcode the exact sequence of moves needed
    moves = [
        # First clear stack 4
        "Move D from 4 to 3",  # Move D (top of stack 4) to stack 3
        "Move B from 4 to 1",  # Move B to join A in stack 1
        "Move I from 4 to 3",  # Move I to stack 3
        "Move J from 4 to 2",  # Move J to stack 2
        "Move H from 4 to 3",  # Move H to stack 3
        # Then arrange other blocks
        "Move G from 1 to 2",  # Move G to stack 2
        "Move E from 2 to 3",  # Move E to stack 3
        "Move F from 2 to 3",  # Move F to stack 3
        "Move C from 3 to 2"   # Move C to stack 2
    ]
    
    # Print the moves
    print('\n'.join(moves))

# Execute the solution
execute_moves()