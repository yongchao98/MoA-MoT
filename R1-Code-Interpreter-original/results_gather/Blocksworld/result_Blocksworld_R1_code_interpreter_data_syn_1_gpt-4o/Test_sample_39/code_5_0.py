# Initial state
stacks = {
    1: ['H', 'B', 'F'],
    2: ['I'],
    3: ['G', 'C', 'D'],
    4: ['J', 'A', 'E']
}

# Corrected moves to achieve the goal state
moves = [
    (1, 3),  # Move F from stack1 to stack3
    (1, 3),  # Move B from stack1 to stack3
    (1, 3),  # Move H from stack1 to stack3
    (3, 1),  # Move D from stack3 to stack1
    (3, 1),  # Move C from stack3 to stack1
    (3, 1),  # Move G from stack3 to stack1
    (4, 2),  # Move E from stack4 to stack2
    (4, 1),  # Move A from stack4 to stack1
    (4, 2),  # Move J from stack4 to stack2
    (2, 3)   # Move I from stack2 to stack3
]

# Execute moves
move_instructions = []
for source, destination in moves:
    block = stacks[source].pop()
    stacks[destination].append(block)
    move_instructions.append(f"Move {block} from {source} to {destination}")

# Print the move instructions
print("\n".join(move_instructions))