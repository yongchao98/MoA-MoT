# Initial state
stacks = {
    1: ['A', 'D'],
    2: [],
    3: ['F', 'B'],
    4: ['E', 'C']
}

# Corrected moves to achieve the goal state
moves = [
    (4, 2),  # Move C from stack4 to stack2
    (4, 1),  # Move E from stack4 to stack1
    (1, 3),  # Move D from stack1 to stack3
    (1, 2),  # Move A from stack1 to stack2
    (3, 2),  # Move F from stack3 to stack2
    (3, 2),  # Move B from stack3 to stack2
    (1, 3)   # Move E from stack1 to stack3
]

# Execute the moves
move_instructions = []
for source, destination in moves:
    block = stacks[source].pop()
    stacks[destination].append(block)
    move_instructions.append(f"Move {block} from {source} to {destination}")

# Print the move instructions
print("\n".join(move_instructions))