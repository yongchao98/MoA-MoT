def print_solution():
    # Initial state
    stacks = [['H', 'G', 'A', 'C'], ['E'], ['F'], ['B', 'D']]
    goal = [['A', 'B', 'E', 'G'], ['C', 'F'], ['D', 'H'], []]
    moves = []
    
    # Helper function to move a block
    def move(from_idx, to_idx):
        if not stacks[from_idx]:
            return False
        block = stacks[from_idx].pop()
        stacks[to_idx].append(block)
        moves.append(f"Move {block} from {from_idx + 1} to {to_idx + 1}")
        return True

    # Step 1: Clear top of stack 1 to get to A
    move(0, 3)  # Move C to stack 4
    move(0, 3)  # Move A to stack 4
    move(0, 2)  # Move G to stack 3
    move(0, 2)  # Move H to stack 3

    # Step 2: Build stack 1 (A B E G)
    move(3, 1)  # Move A back to stack 2
    move(4, 0)  # Move B to stack 1
    move(1, 0)  # Move E to stack 1
    move(2, 0)  # Move G to stack 1

    # Step 3: Build stack 2 (C F)
    move(3, 1)  # Move C to stack 2
    move(2, 1)  # Move F to stack 2

    # Step 4: Build stack 3 (D H)
    move(2, 2)  # Move H to stack 3
    move(3, 2)  # Move D to stack 3

    print('\n'.join(moves))

print_solution()