# Initial state
stacks = {
    1: ['B', 'A'],
    2: ['C', 'F', 'G'],
    3: ['D', 'I', 'E', 'J', 'H']
}

# Moves to achieve the goal state
moves = [
    (3, 2),  # Move H from stack3 to stack2
    (3, 2),  # Move J from stack3 to stack2
    (3, 1),  # Move E from stack3 to stack1
    (3, 1),  # Move I from stack3 to stack1
    (2, 3),  # Move G from stack2 to stack3
    (2, 3),  # Move F from stack2 to stack3
    (2, 3),  # Move C from stack2 to stack3
    (1, 2),  # Move B from stack1 to stack2
    (3, 2),  # Move C from stack3 to stack2
    (2, 2)   # Move J from stack2 to stack2
]

# Execute moves
for move in moves:
    source, destination = move
    block = stacks[source].pop()
    stacks[destination].append(block)
    print(f"Move {block} from {source} to {destination}")