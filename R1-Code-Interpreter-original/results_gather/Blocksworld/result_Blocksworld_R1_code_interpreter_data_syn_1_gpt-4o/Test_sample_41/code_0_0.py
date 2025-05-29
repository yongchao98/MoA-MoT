# Initial state
stacks = {
    1: ['I'],
    2: ['G', 'C', 'H', 'A', 'E'],
    3: ['F', 'B'],
    4: ['D', 'K', 'J']
}

# Moves to achieve the goal state
moves = [
    (2, 3),  # Move E from stack2 to stack3
    (2, 1),  # Move A from stack2 to stack1
    (2, 3),  # Move H from stack2 to stack3
    (2, 3),  # Move C from stack2 to stack3
    (2, 3),  # Move G from stack2 to stack3
    (4, 1),  # Move J from stack4 to stack1
    (4, 3),  # Move K from stack4 to stack3
    (4, 1),  # Move D from stack4 to stack1
    (3, 2),  # Move B from stack3 to stack2
    (3, 1),  # Move F from stack3 to stack1
    (3, 2),  # Move C from stack3 to stack2
    (3, 2),  # Move E from stack3 to stack2
    (3, 2),  # Move G from stack3 to stack2
    (3, 2)   # Move H from stack3 to stack2
]

# Execute the moves
for move in moves:
    source, destination = move
    block = stacks[source].pop()
    stacks[destination].append(block)
    print(f"Move {block} from {source} to {destination}")