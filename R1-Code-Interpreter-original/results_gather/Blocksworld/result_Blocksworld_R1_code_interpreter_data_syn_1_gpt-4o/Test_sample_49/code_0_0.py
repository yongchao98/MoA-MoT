# Initial stacks
stack1 = []
stack2 = ['E', 'G', 'H', 'D', 'A', 'C']
stack3 = []
stack4 = ['B', 'F']

# Moves list
moves = []

# Function to move a block
def move_block(source, destination, source_stack, destination_stack):
    block = source_stack.pop()
    destination_stack.append(block)
    moves.append(f"Move {block} from {source} to {destination}")

# Move C from stack2 to stack1
move_block(2, 1, stack2, stack1)

# Move A from stack2 to stack1
move_block(2, 1, stack2, stack1)

# Move D from stack2 to stack3
move_block(2, 3, stack2, stack3)

# Move H from stack2 to stack1
move_block(2, 1, stack2, stack1)

# Move G from stack2 to stack1
move_block(2, 1, stack2, stack1)

# Move F from stack4 to stack1
move_block(4, 1, stack4, stack1)

# Move B from stack4 to stack3
move_block(4, 3, stack4, stack3)

# Move E from stack2 to stack3
move_block(2, 3, stack2, stack3)

# Print the moves
print("\n".join(moves))