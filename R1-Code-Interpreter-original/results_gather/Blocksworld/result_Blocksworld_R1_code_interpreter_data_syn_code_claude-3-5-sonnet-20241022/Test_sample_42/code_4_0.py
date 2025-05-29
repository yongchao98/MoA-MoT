def make_move(stacks, from_stack, to_stack):
    if not stacks[from_stack]:
        return None
    block = stacks[from_stack].pop()
    stacks[to_stack].append(block)
    return f"Move {block} from {from_stack + 1} to {to_stack + 1}"

def print_state(stacks):
    print(f"Stack1: {stacks[0]}")
    print(f"Stack2: {stacks[1]}")
    print(f"Stack3: {stacks[2]}\n")

def solve_blocksworld():
    # Initial state
    stacks = [['A'], [], ['B', 'C', 'D', 'E']]
    moves = []
    
    # Step 1: First get E out of the way
    moves.append(make_move(stacks, 2, 1))  # Move E from 3 to 2
    
    # Step 2: Get D to its final position
    moves.append(make_move(stacks, 2, 0))  # Move D from 3 to 1
    
    # Step 3: Move A to stack3 (it needs to be at bottom)
    moves.append(make_move(stacks, 0, 2))  # Move A from 1 to 3
    
    # Step 4: Move C to its final position
    moves.append(make_move(stacks, 2, 2))  # Move C from 3 to 3
    
    # Step 5: Move B to stack2 with E
    moves.append(make_move(stacks, 2, 1))  # Move B from 3 to 2
    
    print("<<<" + "\n".join([move for move in moves if move is not None]) + ">>>")

solve_blocksworld()