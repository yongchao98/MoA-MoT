def is_valid_move(stacks, source, dest):
    # Check if source stack is not empty
    if not stacks[source-1]:
        return False
    return True

def make_move(stacks, source, dest):
    if not is_valid_move(stacks, source, dest):
        return False
    block = stacks[source-1].pop()
    stacks[dest-1].append(block)
    return True

def print_state(stacks):
    for i, stack in enumerate(stacks):
        print(f"Stack {i+1}: {' '.join(stack)}")

def find_solution():
    # Initial state
    stacks = [
        ['D'],           # stack 1
        ['F', 'A', 'B'], # stack 2
        ['E'],          # stack 3
        ['C']           # stack 4
    ]
    
    # Goal state
    goal = [
        [],                 # stack 1
        ['B', 'C', 'D', 'E'], # stack 2
        ['A', 'F']         # stack 3
    ]
    
    moves = []
    
    # Move B from stack2 to stack1 (temporary)
    moves.append(('B', 2, 1))
    make_move(stacks, 2, 1)
    
    # Move A from stack2 to stack4 (temporary)
    moves.append(('A', 2, 4))
    make_move(stacks, 2, 4)
    
    # Move F from stack2 to stack3
    moves.append(('F', 2, 3))
    make_move(stacks, 2, 3)
    
    # Move A from stack4 to stack3
    moves.append(('A', 4, 3))
    make_move(stacks, 4, 3)
    
    # Move C from stack4 to stack2
    moves.append(('C', 4, 2))
    make_move(stacks, 4, 2)
    
    # Move B from stack1 to stack2
    moves.append(('B', 1, 2))
    make_move(stacks, 1, 2)
    
    # Move D from stack1 to stack2
    moves.append(('D', 1, 2))
    make_move(stacks, 1, 2)
    
    # Move E from stack3 to stack2
    moves.append(('E', 3, 2))
    make_move(stacks, 3, 2)
    
    # Format and print the moves
    formatted_moves = []
    for block, source, dest in moves:
        formatted_moves.append(f"Move {block} from {source} to {dest}")
    
    print("<<<" + "\n".join(formatted_moves) + ">>>")

find_solution()