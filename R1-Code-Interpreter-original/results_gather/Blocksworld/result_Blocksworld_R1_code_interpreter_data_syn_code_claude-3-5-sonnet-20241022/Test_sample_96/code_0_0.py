def is_valid_move(stacks, from_stack, to_stack):
    # Check if source stack is not empty
    if not stacks[from_stack-1]:
        return False
    return True

def make_move(stacks, from_stack, to_stack):
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    return block

def print_state(stacks):
    for i, stack in enumerate(stacks, 1):
        print(f"Stack {i}: {' '.join(stack)}")
    print()

def find_solution():
    # Initial state
    stacks = [
        ['E', 'B', 'J', 'I'],
        ['A', 'F', 'K', 'C'],
        ['D', 'G', 'H']
    ]
    
    # Goal state
    goal = [
        ['A', 'C', 'F', 'J', 'K'],
        ['D', 'G'],
        ['B', 'E', 'H', 'I']
    ]
    
    moves = []
    
    # Manual solution (based on analysis)
    # Move I from stack1 to stack3
    moves.append(('I', 1, 3))
    stacks[0].pop()
    stacks[2].append('I')
    
    # Move J from stack1 to stack1 (temporary)
    moves.append(('J', 1, 1))
    
    # Move C from stack2 to stack1
    moves.append(('C', 2, 1))
    stacks[1].pop()
    stacks[0].append('C')
    
    # Move K from stack2 to stack1
    moves.append(('K', 2, 1))
    stacks[1].pop()
    stacks[0].append('K')
    
    # Move F from stack2 to stack1
    moves.append(('F', 2, 1))
    stacks[1].pop()
    stacks[0].append('F')
    
    # Move A from stack2 to stack1
    moves.append(('A', 2, 1))
    stacks[1].pop()
    stacks[0].append('A')
    
    # Move H from stack3 to stack2
    moves.append(('H', 3, 2))
    stacks[2].pop()
    stacks[1].append('H')
    
    # Move G from stack3 to stack2
    moves.append(('G', 3, 2))
    stacks[2].pop()
    stacks[1].append('G')
    
    # Move D from stack3 to stack2
    moves.append(('D', 3, 2))
    stacks[2].pop()
    stacks[1].append('D')
    
    # Move B from stack1 to stack3
    moves.append(('B', 1, 3))
    stacks[0].remove('B')
    stacks[2].append('B')
    
    # Move E from stack1 to stack3
    moves.append(('E', 1, 3))
    stacks[0].remove('E')
    stacks[2].append('E')
    
    # Format and print the moves
    for block, source, dest in moves:
        print(f"Move {block} from {source} to {dest}")

find_solution()