# Initial state
stacks = {
    1: ['J'],
    2: ['C', 'K', 'F', 'E', 'B', 'I'],
    3: ['G', 'A'],
    4: ['H', 'D']
}

# Goal state
goal = {
    1: ['F', 'G', 'H', 'I', 'J'],
    2: ['C', 'D', 'E', 'K'],
    3: ['A', 'B']
}

# Function to move a block
def move_block(stacks, source, destination):
    block = stacks[source].pop()
    stacks[destination].append(block)
    return f"Move {block} from {source} to {destination}"

# List to store the sequence of moves
moves = []

# Move sequence to achieve the goal state
moves.append(move_block(stacks, 2, 4))  # Move I from 2 to 4
moves.append(move_block(stacks, 2, 1))  # Move B from 2 to 1
moves.append(move_block(stacks, 2, 3))  # Move E from 2 to 3
moves.append(move_block(stacks, 2, 3))  # Move F from 2 to 3
moves.append(move_block(stacks, 2, 4))  # Move K from 2 to 4
moves.append(move_block(stacks, 2, 3))  # Move C from 2 to 3
moves.append(move_block(stacks, 4, 2))  # Move K from 4 to 2
moves.append(move_block(stacks, 4, 1))  # Move I from 4 to 1
moves.append(move_block(stacks, 3, 2))  # Move C from 3 to 2
moves.append(move_block(stacks, 3, 1))  # Move F from 3 to 1
moves.append(move_block(stacks, 3, 1))  # Move E from 3 to 1
moves.append(move_block(stacks, 4, 3))  # Move D from 4 to 3
moves.append(move_block(stacks, 4, 1))  # Move H from 4 to 1
moves.append(move_block(stacks, 3, 1))  # Move G from 3 to 1
moves.append(move_block(stacks, 2, 3))  # Move B from 2 to 3
moves.append(move_block(stacks, 2, 3))  # Move C from 2 to 3

# Print the sequence of moves
for move in moves:
    print(move)