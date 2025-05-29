# Initial state
stacks = {
    1: [],
    2: ['E', 'B', 'I'],
    3: ['F', 'D'],
    4: ['C', 'H', 'G', 'A']
}

# Function to move a block
def move_block(stacks, source, destination):
    block = stacks[source].pop()
    stacks[destination].append(block)
    return f"Move {block} from {source} to {destination}"

# List to store the moves
moves = []

# Move sequence to achieve the goal state
moves.append(move_block(stacks, 4, 1))  # Move A from 4 to 1
moves.append(move_block(stacks, 4, 2))  # Move G from 4 to 2
moves.append(move_block(stacks, 4, 3))  # Move H from 4 to 3
moves.append(move_block(stacks, 4, 3))  # Move C from 4 to 3
moves.append(move_block(stacks, 2, 3))  # Move I from 2 to 3
moves.append(move_block(stacks, 2, 1))  # Move B from 2 to 1
moves.append(move_block(stacks, 3, 1))  # Move D from 3 to 1
moves.append(move_block(stacks, 3, 2))  # Move C from 3 to 2
moves.append(move_block(stacks, 3, 2))  # Move H from 3 to 2
moves.append(move_block(stacks, 2, 3))  # Move G from 2 to 3
moves.append(move_block(stacks, 2, 3))  # Move E from 2 to 3
moves.append(move_block(stacks, 2, 3))  # Move F from 2 to 3

# Print the moves
for move in moves:
    print(move)