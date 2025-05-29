# Initial configuration
stacks = {
    1: ['C'],
    2: ['B', 'E', 'A', 'J', 'F'],
    3: ['G', 'I'],
    4: ['D', 'H']
}

# Function to move a block
def move_block(stacks, source, destination):
    block = stacks[source].pop()
    stacks[destination].append(block)
    return f"Move {block} from {source} to {destination}"

# List to store the moves
moves = []

# Move sequence to achieve the goal state
moves.append(move_block(stacks, 2, 3))  # Move F from 2 to 3
moves.append(move_block(stacks, 2, 1))  # Move J from 2 to 1
moves.append(move_block(stacks, 2, 3))  # Move A from 2 to 3
moves.append(move_block(stacks, 2, 4))  # Move E from 2 to 4
moves.append(move_block(stacks, 2, 4))  # Move B from 2 to 4
moves.append(move_block(stacks, 1, 2))  # Move C from 1 to 2
moves.append(move_block(stacks, 4, 2))  # Move B from 4 to 2
moves.append(move_block(stacks, 4, 2))  # Move E from 4 to 2
moves.append(move_block(stacks, 4, 2))  # Move H from 4 to 2
moves.append(move_block(stacks, 4, 2))  # Move D from 4 to 2
moves.append(move_block(stacks, 3, 2))  # Move A from 3 to 2
moves.append(move_block(stacks, 3, 1))  # Move I from 3 to 1
moves.append(move_block(stacks, 3, 2))  # Move G from 3 to 2
moves.append(move_block(stacks, 3, 2))  # Move F from 3 to 2
moves.append(move_block(stacks, 1, 3))  # Move J from 1 to 3

# Print the moves
print("\n".join(moves))