# Initial state
stacks = {
    1: ['G', 'I', 'D'],
    2: ['E', 'A', 'K'],
    3: ['C', 'B', 'H', 'J', 'F']
}

# Corrected moves to achieve the goal state
moves = [
    (3, 1),  # Move F from stack3 to stack1
    (3, 2),  # Move J from stack3 to stack2
    (3, 2),  # Move H from stack3 to stack2
    (3, 2),  # Move B from stack3 to stack2
    (2, 1),  # Move A from stack2 to stack1
    (2, 3),  # Move E from stack2 to stack3
    (1, 2),  # Move I from stack1 to stack2
    (1, 3)   # Move G from stack1 to stack3
]

# Execute moves
move_instructions = []
for source, destination in moves:
    block = stacks[source].pop()
    stacks[destination].append(block)
    move_instructions.append(f"Move {block} from {source} to {destination}")

# Print the move instructions
print("\n".join(move_instructions))