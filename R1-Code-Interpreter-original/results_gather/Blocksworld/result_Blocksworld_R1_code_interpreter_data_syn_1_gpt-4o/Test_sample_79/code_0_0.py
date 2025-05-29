# Initial State
stacks = {
    1: ['F', 'H'],
    2: ['G', 'D'],
    3: ['E', 'C', 'B'],
    4: ['A', 'I']
}

# Moves to reach the goal state
moves = []

# Move D from stack2 to stack1
stacks[1].append(stacks[2].pop())
moves.append("Move D from 2 to 1")

# Move B from stack3 to stack2
stacks[2].append(stacks[3].pop())
moves.append("Move B from 3 to 2")

# Move C from stack3 to stack2
stacks[2].append(stacks[3].pop())
moves.append("Move C from 3 to 2")

# Move E from stack3 to stack2
stacks[2].append(stacks[3].pop())
moves.append("Move E from 3 to 2")

# Move G from stack2 to stack3
stacks[3].append(stacks[2].pop())
moves.append("Move G from 2 to 3")

# Move H from stack1 to stack3
stacks[3].append(stacks[1].pop())
moves.append("Move H from 1 to 3")

# Move I from stack4 to stack3
stacks[3].append(stacks[4].pop())
moves.append("Move I from 4 to 3")

# Move A from stack4 to stack3
stacks[3].append(stacks[4].pop())
moves.append("Move A from 4 to 3")

# Print the moves
print("\n".join(moves))