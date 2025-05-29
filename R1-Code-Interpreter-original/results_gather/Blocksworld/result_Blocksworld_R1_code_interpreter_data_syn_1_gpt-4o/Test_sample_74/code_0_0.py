# Initial state
stacks = {
    1: ['B', 'G', 'F'],
    2: ['A', 'E', 'C'],
    3: ['D']
}

# Moves to achieve the goal state
moves = []

# Move D from stack3 to stack1
stacks[1].append(stacks[3].pop())
moves.append("Move D from 3 to 1")

# Move F from stack1 to stack3
stacks[3].append(stacks[1].pop())
moves.append("Move F from 1 to 3")

# Move G from stack1 to stack1
stacks[1].append(stacks[1].pop())
moves.append("Move G from 1 to 1")

# Move F from stack3 to stack1
stacks[1].append(stacks[3].pop())
moves.append("Move F from 3 to 1")

# Move C from stack2 to stack3
stacks[3].append(stacks[2].pop())
moves.append("Move C from 2 to 3")

# Move E from stack2 to stack2
stacks[2].append(stacks[2].pop())
moves.append("Move E from 2 to 2")

# Move A from stack2 to stack3
stacks[3].append(stacks[2].pop())
moves.append("Move A from 2 to 3")

# Move B from stack1 to stack3
stacks[3].append(stacks[1].pop())
moves.append("Move B from 1 to 3")

# Print the moves
print("\n".join(moves))