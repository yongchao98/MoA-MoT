# Initial state
stacks = {
    1: ['J', 'B'],
    2: ['F', 'I', 'E'],
    3: ['D', 'A', 'C', 'H', 'G']
}

# Goal state
goal = {
    1: ['B', 'C', 'G', 'H', 'I', 'J'],
    2: ['D', 'F'],
    3: ['A', 'E']
}

# Function to move a block
def move_block(source, destination):
    block = stacks[source].pop()
    stacks[destination].append(block)
    return f"Move {block} from {source} to {destination}"

# List to store the sequence of moves
moves = []

# Move sequence to achieve the goal state
moves.append(move_block(1, 3))  # Move B from 1 to 3
moves.append(move_block(1, 2))  # Move J from 1 to 2
moves.append(move_block(2, 1))  # Move E from 2 to 1
moves.append(move_block(2, 3))  # Move I from 2 to 3
moves.append(move_block(2, 1))  # Move F from 2 to 1
moves.append(move_block(3, 2))  # Move G from 3 to 2
moves.append(move_block(3, 1))  # Move H from 3 to 1
moves.append(move_block(3, 1))  # Move C from 3 to 1
moves.append(move_block(3, 2))  # Move A from 3 to 2
moves.append(move_block(3, 1))  # Move D from 3 to 1
moves.append(move_block(2, 3))  # Move A from 2 to 3
moves.append(move_block(2, 3))  # Move D from 2 to 3
moves.append(move_block(1, 2))  # Move F from 1 to 2
moves.append(move_block(1, 3))  # Move C from 1 to 3
moves.append(move_block(1, 3))  # Move H from 1 to 3
moves.append(move_block(1, 3))  # Move G from 1 to 3
moves.append(move_block(2, 1))  # Move F from 2 to 1
moves.append(move_block(2, 1))  # Move D from 2 to 1
moves.append(move_block(3, 1))  # Move G from 3 to 1
moves.append(move_block(3, 1))  # Move H from 3 to 1
moves.append(move_block(3, 1))  # Move C from 3 to 1
moves.append(move_block(3, 1))  # Move B from 3 to 1
moves.append(move_block(2, 3))  # Move D from 2 to 3
moves.append(move_block(2, 3))  # Move F from 2 to 3

# Print the sequence of moves
print("\n".join(moves))