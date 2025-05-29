# Initial stacks
stacks = {
    1: ['B', 'I'],
    2: ['G', 'F'],
    3: ['C', 'D', 'E'],
    4: ['A', 'H']
}

# Moves to achieve the goal state
moves = [
    (3, 1),  # Move E from stack3 to stack1
    (3, 2),  # Move D from stack3 to stack2
    (2, 3),  # Move F from stack2 to stack3
    (2, 3),  # Move G from stack2 to stack3
    (4, 3),  # Move H from stack4 to stack3
    (4, 2),  # Move A from stack4 to stack2
    (1, 3),  # Move I from stack1 to stack3
    (1, 3)   # Move B from stack1 to stack3
]

# Execute the moves
move_instructions = []
for source, destination in moves:
    block = stacks[source].pop()
    stacks[destination].append(block)
    move_instructions.append(f"Move {block} from {source} to {destination}")

# Print the move instructions
print("\n".join(move_instructions))