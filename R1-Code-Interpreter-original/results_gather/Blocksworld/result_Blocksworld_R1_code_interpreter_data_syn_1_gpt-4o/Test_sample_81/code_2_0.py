# Initial state
stacks = {
    1: ['D'],
    2: ['A', 'J'],
    3: ['H', 'F', 'C', 'E', 'I'],
    4: ['B', 'G']
}

# Corrected moves to achieve the goal state
moves = [
    (1, 2),  # Move D from stack1 to stack2
    (2, 3),  # Move J from stack2 to stack3
    (3, 1),  # Move I from stack3 to stack1
    (3, 3),  # Move E from stack3 to stack3 (to free up F)
    (3, 1),  # Move F from stack3 to stack1
    (3, 2),  # Move C from stack3 to stack2
    (3, 3),  # Move H from stack3 to stack3 (to free up J)
    (3, 3),  # Move J from stack3 to stack3 (to complete stack3)
    (4, 2),  # Move G from stack4 to stack2
    (4, 2)   # Move B from stack4 to stack2
]

# Execute the corrected moves
move_instructions = []
for source, destination in moves:
    block = stacks[source].pop()
    stacks[destination].append(block)
    move_instructions.append(f"Move {block} from {source} to {destination}")

# Print the move instructions
print("\n".join(move_instructions))