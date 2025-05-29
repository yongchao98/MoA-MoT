# Initial stacks
stack1 = ['J']
stack2 = ['D', 'H', 'F', 'I', 'B', 'A', 'G']
stack3 = ['E', 'C']

# Moves list
moves = []

# Function to move a block
def move_block(source, destination, source_stack, destination_stack):
    block = source_stack.pop()
    destination_stack.append(block)
    moves.append(f"Move {block} from {source} to {destination}")

# Move blocks to free G and I
move_block(2, 3, stack2, stack3)  # Move G from stack2 to stack3
move_block(2, 3, stack2, stack3)  # Move A from stack2 to stack3
move_block(2, 3, stack2, stack3)  # Move B from stack2 to stack3
move_block(2, 3, stack2, stack3)  # Move I from stack2 to stack3

# Move G and I to stack1
move_block(3, 1, stack3, stack1)  # Move I from stack3 to stack1
move_block(3, 1, stack3, stack1)  # Move G from stack3 to stack1

# Move blocks to achieve the goal state
move_block(3, 2, stack3, stack2)  # Move B from stack3 to stack2
move_block(3, 2, stack3, stack2)  # Move A from stack3 to stack2
move_block(3, 2, stack3, stack2)  # Move G from stack3 to stack2
move_block(3, 2, stack3, stack2)  # Move I from stack3 to stack2
move_block(1, 3, stack1, stack3)  # Move J from stack1 to stack3

# Print the moves
print("\n".join(moves))