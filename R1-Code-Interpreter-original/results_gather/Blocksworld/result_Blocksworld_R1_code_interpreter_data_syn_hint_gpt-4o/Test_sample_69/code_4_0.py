# Initial stacks
stacks = {
    "stack1": ["H", "I", "B", "E", "A"],
    "stack2": ["J"],
    "stack3": ["G"],
    "stack4": ["F", "D", "K", "C"]
}

# Function to move a block from one stack to another
def move_block(source, destination):
    block = stacks[source].pop()
    stacks[destination].append(block)
    return f"Move {block} from {source} to {destination}"

# List to store the sequence of moves
moves = []

# Perform the moves
moves.append(move_block("stack1", "stack3"))  # Move A
moves.append(move_block("stack1", "stack3"))  # Move E
moves.append(move_block("stack1", "stack2"))  # Move B
moves.append(move_block("stack1", "stack2"))  # Move I
moves.append(move_block("stack4", "stack2"))  # Move C
moves.append(move_block("stack4", "stack3"))  # Move K
moves.append(move_block("stack4", "stack3"))  # Move D
moves.append(move_block("stack4", "stack2"))  # Move F
moves.append(move_block("stack3", "stack2"))  # Move G
moves.append(move_block("stack2", "stack3"))  # Move J

# Print the sequence of moves
for move in moves:
    print(move)