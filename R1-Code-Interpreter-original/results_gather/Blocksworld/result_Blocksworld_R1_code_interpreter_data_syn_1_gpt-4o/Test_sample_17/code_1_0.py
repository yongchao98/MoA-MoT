# Initial state
stacks = {
    'stack1': ['C', 'F', 'E'],
    'stack2': ['D'],
    'stack3': [],
    'stack4': ['A', 'B']
}

# Corrected moves to achieve the goal state
moves = [
    ('B', 'stack4', 'stack2'),
    ('A', 'stack4', 'stack3'),
    ('E', 'stack1', 'stack3'),
    ('F', 'stack1', 'stack4'),
    ('E', 'stack3', 'stack1'),
    ('F', 'stack4', 'stack1')
]

# Function to perform a move
def move_block(block, source, destination):
    if stacks[source][-1] == block:
        stacks[source].pop()
        stacks[destination].append(block)
    else:
        raise ValueError(f"Cannot move block {block} from {source} as it is not on top.")

# Execute the moves
for block, source, destination in moves:
    move_block(block, source, destination)

# Print the final state of the stacks
print(stacks)