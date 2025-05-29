def print_stacks(stacks):
    for i, s in enumerate(stacks):
        print(f"Stack {i+1}: {' '.join(s)}")

def get_top(stack):
    return stack[-1] if stack else None

def make_move(stacks, from_stack, to_stack, moves):
    # Convert to 0-based indexing
    from_idx = from_stack - 1
    to_idx = to_stack - 1
    
    # Verify source stack has blocks
    if not stacks[from_idx]:
        return False
    
    # Get the top block without removing it first
    block = stacks[from_idx][-1]
    
    # Now actually make the move
    stacks[from_idx].pop()
    stacks[to_idx].append(block)
    
    # Record the move
    moves.append(f"Move {block} from {from_stack} to {to_stack}")
    return True

def solve_blocks():
    # Initialize stacks
    stacks = [
        ['F'],           # stack 1
        ['B'],           # stack 2
        ['D', 'A', 'E', 'C']  # stack 3
    ]
    
    moves = []
    
    # Explicitly planned moves
    # Move C from stack 3 to stack 2
    stacks[2].pop()  # Remove C
    stacks[1].append('C')
    moves.append("Move C from 3 to 2")
    
    # Move E from stack 3 to stack 1
    stacks[2].pop()  # Remove E
    stacks[0].append('E')
    moves.append("Move E from 3 to 1")
    
    # Move A from stack 3 to stack 1
    stacks[2].pop()  # Remove A
    stacks[0].append('A')
    moves.append("Move A from 3 to 1")
    
    # Move B from stack 2 to stack 3
    stacks[1].pop(0)  # Remove B
    stacks[2].append('B')
    moves.append("Move B from 2 to 3")
    
    # Move A from stack 1 to stack 3
    stacks[0].pop()  # Remove A
    stacks[2].append('A')
    moves.append("Move A from 1 to 3")
    
    # Move E from stack 1 to stack 2
    stacks[0].pop()  # Remove E
    stacks[1].append('E')
    moves.append("Move E from 1 to 2")
    
    print("<<<" + "\n".join(moves) + ">>>")

solve_blocks()