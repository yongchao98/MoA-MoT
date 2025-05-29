# Initial state
stack1 = ['D', 'E', 'H', 'G']
stack2 = ['A', 'F', 'I']
stack3 = ['B', 'C']

# Moves to achieve the goal state
moves = []

# Move I from stack2 to stack1
stack1.append(stack2.pop())
moves.append("Move I from 2 to 1")

# Move H from stack1 to stack2
stack2.append(stack1.pop())
moves.append("Move H from 1 to 2")

# Move E from stack1 to stack3
stack3.append(stack1.pop())
moves.append("Move E from 1 to 3")

# Move D from stack1 to stack2
stack2.append(stack1.pop())
moves.append("Move D from 1 to 2")

# Move C from stack3 to stack2
stack2.append(stack3.pop())
moves.append("Move C from 3 to 2")

# Print the moves
for move in moves:
    print(move)