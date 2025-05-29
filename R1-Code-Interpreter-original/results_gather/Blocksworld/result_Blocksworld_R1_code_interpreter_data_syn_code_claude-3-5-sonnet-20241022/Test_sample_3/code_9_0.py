def print_state(stacks):
    for i, stack in enumerate(stacks, 1):
        print(f"stack{i}: {' '.join(stack)}")
    print()

def is_valid_move(stacks, from_stack, to_stack):
    if not stacks[from_stack-1]:
        return False
    return True

def make_move(stacks, from_stack, to_stack):
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    return block

def find_solution():
    # Initial state
    stacks = [['A', 'D'], [], ['F', 'B'], ['E', 'C']]
    moves = []
    
    print("Initial state:")
    print_state(stacks)
    
    # Move C from stack4 to stack2
    if is_valid_move(stacks, 4, 2):
        block = make_move(stacks, 4, 2)
        moves.append(f"Move {block} from 4 to 2")
        print_state(stacks)
    
    # Move B from stack3 to stack2
    if is_valid_move(stacks, 3, 2):
        block = make_move(stacks, 3, 2)
        moves.append(f"Move {block} from 3 to 2")
        print_state(stacks)
    
    # Move D from stack1 to stack3
    if is_valid_move(stacks, 1, 3):
        block = make_move(stacks, 1, 3)
        moves.append(f"Move {block} from 1 to 3")
        print_state(stacks)
    
    # Move E from stack4 to stack1
    if is_valid_move(stacks, 4, 1):
        block = make_move(stacks, 4, 1)
        moves.append(f"Move {block} from 4 to 1")
        print_state(stacks)
    
    # Move F from stack3 to stack4 (temporary)
    if is_valid_move(stacks, 3, 4):
        block = make_move(stacks, 3, 4)
        moves.append(f"Move {block} from 3 to 4")
        print_state(stacks)
    
    # Move A from stack1 to stack2
    if is_valid_move(stacks, 1, 2):
        block = make_move(stacks, 1, 2)
        moves.append(f"Move {block} from 1 to 2")
        print_state(stacks)
    
    # Move F from stack4 to stack2
    if is_valid_move(stacks, 4, 2):
        block = make_move(stacks, 4, 2)
        moves.append(f"Move {block} from 4 to 2")
        print_state(stacks)
    
    print("<<<" + "\n".join(moves) + ">>>")

find_solution()