# Initial configuration
stacks = {
    1: ['I', 'D', 'B'],
    2: ['G', 'J', 'A', 'C', 'E'],
    3: ['H', 'F']
}

# Goal configuration
goal = {
    1: ['E', 'H'],
    2: ['A', 'G'],
    3: ['B', 'C', 'D', 'F', 'I', 'J']
}

# Function to move a block
def move_block(stacks, source, destination):
    block = stacks[source].pop()
    stacks[destination].append(block)
    return f"Move {block} from {source} to {destination}"

# List to store the moves
moves = []

# Move sequence to achieve the goal state
moves.append(move_block(stacks, 1, 3))  # Move B from 1 to 3
moves.append(move_block(stacks, 1, 3))  # Move D from 1 to 3
moves.append(move_block(stacks, 1, 3))  # Move I from 1 to 3
moves.append(move_block(stacks, 2, 1))  # Move E from 2 to 1
moves.append(move_block(stacks, 3, 1))  # Move F from 3 to 1
moves.append(move_block(stacks, 3, 2))  # Move I from 3 to 2
moves.append(move_block(stacks, 3, 2))  # Move D from 3 to 2
moves.append(move_block(stacks, 3, 2))  # Move B from 3 to 2
moves.append(move_block(stacks, 2, 3))  # Move C from 2 to 3
moves.append(move_block(stacks, 2, 3))  # Move A from 2 to 3
moves.append(move_block(stacks, 2, 3))  # Move J from 2 to 3
moves.append(move_block(stacks, 2, 3))  # Move G from 2 to 3
moves.append(move_block(stacks, 1, 3))  # Move F from 1 to 3
moves.append(move_block(stacks, 1, 3))  # Move E from 1 to 3
moves.append(move_block(stacks, 3, 1))  # Move H from 3 to 1
moves.append(move_block(stacks, 3, 2))  # Move E from 3 to 2

# Print the moves
for move in moves:
    print(move)