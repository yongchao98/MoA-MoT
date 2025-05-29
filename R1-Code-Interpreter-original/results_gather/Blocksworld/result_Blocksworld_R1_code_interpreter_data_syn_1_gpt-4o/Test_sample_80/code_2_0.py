# Initial state
stack1 = ['E', 'C', 'J']
stack2 = ['H', 'D']
stack3 = ['G', 'I', 'A', 'B', 'F']

# Moves to achieve the goal state
moves = [
    ('F', 3, 1),
    ('B', 3, 2),
    ('A', 3, 1),
    ('G', 3, 1),
    ('D', 2, 1),
    ('J', 1, 2)
]

# Function to perform a move
def move_block(source, destination):
    block = source.pop()
    destination.append(block)

# Perform the moves
for block, src, dest in moves:
    if src == 1:
        source = stack1
    elif src == 2:
        source = stack2
    else:
        source = stack3

    if dest == 1:
        destination = stack1
    elif dest == 2:
        destination = stack2
    else:
        destination = stack3

    move_block(source, destination)

# Print the final stacks
print("stack1:", stack1)
print("stack2:", stack2)
print("stack3:", stack3)