def print_solution():
    # Initial state
    stacks = [['H', 'G', 'A', 'C'], ['E'], ['F'], ['B', 'D']]
    moves = []
    
    # Helper function to move a block
    def move(from_idx, to_idx):
        if 0 <= from_idx < len(stacks) and 0 <= to_idx < len(stacks) and stacks[from_idx]:
            block = stacks[from_idx].pop()
            stacks[to_idx].append(block)
            moves.append(f"Move {block} from {from_idx + 1} to {to_idx + 1}")
            return True
        return False

    # First clear C from top of stack 1
    move(0, 1)  # Move C from stack 1 to 2
    
    # Move A temporarily
    move(0, 2)  # Move A from stack 1 to 3
    
    # Move G temporarily
    move(0, 3)  # Move G from stack 1 to 4
    
    # Move H to its final position
    move(0, 2)  # Move H from stack 1 to 3
    
    # Move D to its final position
    move(3, 2)  # Move D from stack 4 to 3
    
    # Start building stack 1
    move(3, 0)  # Move B from stack 4 to 1
    move(1, 0)  # Move E from stack 2 to 1
    move(3, 0)  # Move G from stack 4 to 1
    
    # Build stack 2
    move(1, 1)  # Move C to stack 2
    move(2, 1)  # Move F to stack 2

    print('\n'.join(moves))

print_solution()