# Initial stacks
stacks = {
    1: ['C', 'B'],
    2: ['D', 'H', 'E'],
    3: ['F', 'A'],
    4: ['G']
}

# Function to move a block from one stack to another
def move_block(source, destination):
    block = stacks[source].pop()
    stacks[destination].append(block)
    return f"Move {block} from {source} to {destination}"

# List to store the moves
moves = []

# Execute the planned moves
moves.append(move_block(2, 1))  # Move E from stack2 to stack1
moves.append(move_block(2, 1))  # Move H from stack2 to stack1
moves.append(move_block(3, 1))  # Move A from stack3 to stack1
moves.append(move_block(4, 3))  # Move G from stack4 to stack3
moves.append(move_block(1, 2))  # Move C from stack1 to stack2
moves.append(move_block(1, 2))  # Move B from stack1 to stack2
moves.append(move_block(2, 1))  # Move D from stack2 to stack1
moves.append(move_block(2, 1))  # Move B from stack2 to stack1
moves.append(move_block(2, 1))  # Move C from stack2 to stack1

# Print the moves
for move in moves:
    print(move)