# Initial state
stacks = {
    1: ['G', 'C'],
    2: ['D', 'A', 'E', 'B'],
    3: ['I', 'H', 'F']
}

# Moves to achieve the goal state
moves = []

# Move B from stack2 to stack1
stacks[1].append(stacks[2].pop())
moves.append("Move B from 2 to 1")

# Move A from stack2 to stack3
stacks[3].append(stacks[2].pop())
moves.append("Move A from 2 to 3")

# Move D from stack2 to stack3
stacks[3].append(stacks[2].pop())
moves.append("Move D from 2 to 3")

# Move C from stack1 to stack3
stacks[3].append(stacks[1].pop())
moves.append("Move C from 1 to 3")

# Move G from stack1 to stack3
stacks[3].append(stacks[1].pop())
moves.append("Move G from 1 to 3")

# Move I from stack3 to stack1
stacks[1].append(stacks[3].pop())
moves.append("Move I from 3 to 1")

# Print the moves
print("\n".join(moves))