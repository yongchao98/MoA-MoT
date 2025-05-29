# Initial state
stacks = {
    1: ['D'],
    2: ['E', 'A'],
    3: ['C', 'B'],
    4: []
}

# Function to move a block from one stack to another
def move_block(source, destination):
    block = stacks[source].pop()
    stacks[destination].append(block)
    return f"Move {block} from {source} to {destination}"

# List to store the sequence of moves
moves = []

# Execute the plan
moves.append(move_block(1, 4))  # Move D from stack1 to stack4
moves.append(move_block(2, 1))  # Move A from stack2 to stack1
moves.append(move_block(2, 3))  # Move E from stack2 to stack3
moves.append(move_block(3, 2))  # Move B from stack3 to stack2
moves.append(move_block(3, 1))  # Move C from stack3 to stack1
moves.append(move_block(3, 2))  # Move E from stack3 to stack2
moves.append(move_block(4, 2))  # Move D from stack4 to stack2

# Print the sequence of moves
print("\n".join(moves))