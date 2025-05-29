# Initial state
stacks = {
    1: ['E'],
    2: ['B'],
    3: ['H', 'C', 'A', 'F'],
    4: ['D', 'G']
}

# Moves to achieve the goal state
moves = []

# Move F from stack3 to stack2
stacks[2].append(stacks[3].pop())
moves.append("Move F from 3 to 2")

# Move A from stack3 to stack3 (no actual move needed, just for clarity)
# stacks[3].append(stacks[3].pop())
moves.append("Move A from 3 to 3")

# Move G from stack4 to stack1
stacks[1].append(stacks[4].pop())
moves.append("Move G from 4 to 1")

# Move D from stack4 to stack1
stacks[1].append(stacks[4].pop())
moves.append("Move D from 4 to 1")

# Move C from stack3 to stack2
stacks[2].append(stacks[3].pop())
moves.append("Move C from 3 to 2")

# Move H from stack3 to stack1
stacks[1].append(stacks[3].pop())
moves.append("Move H from 3 to 1")

# Print the moves
print("\n".join(moves))