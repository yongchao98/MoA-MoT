def print_state(stacks):
    for i, stack in enumerate(stacks):
        print(f"Stack {i+1}: {stack}")

def make_move(stacks, from_stack, to_stack, moves):
    if not stacks[from_stack-1]:  # Check if source stack is empty
        return False
    block = stacks[from_stack-1].pop()
    stacks[to_stack-1].append(block)
    moves.append(f"Move {block} from {from_stack} to {to_stack}")
    return True

def find_solution():
    # Initial state
    stacks = [
        ['F'],           # stack 1
        ['B'],           # stack 2
        ['D', 'A', 'E', 'C']  # stack 3
    ]
    
    moves = []
    
    # Step 1: Move C from stack 3 to stack 2
    make_move(stacks, 3, 2, moves)
    
    # Step 2: Move E from stack 3 to stack 1
    make_move(stacks, 3, 1, moves)
    
    # Step 3: Move A from stack 3 to stack 1
    make_move(stacks, 3, 1, moves)
    
    # Step 4: Move B from stack 2 to stack 3
    make_move(stacks, 2, 3, moves)
    
    # Step 5: Move A from stack 1 to stack 3
    make_move(stacks, 1, 3, moves)
    
    # Step 6: Move E from stack 1 to stack 2
    make_move(stacks, 1, 2, moves)
    
    print("<<<" + "\n".join(moves) + ">>>")

find_solution()