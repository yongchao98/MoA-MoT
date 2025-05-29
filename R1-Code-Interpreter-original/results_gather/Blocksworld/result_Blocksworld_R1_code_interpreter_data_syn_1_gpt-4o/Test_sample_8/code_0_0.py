# Initial state
stack1 = ['D']
stack2 = ['G', 'C', 'A', 'F']
stack3 = ['B', 'I', 'E', 'H']

# Moves to achieve the goal state
moves = []

# Move H from stack3 to stack1
stack1.append(stack3.pop())
moves.append("Move H from 3 to 1")

# Move E from stack3 to stack2
stack2.append(stack3.pop())
moves.append("Move E from 3 to 2")

# Move I from stack3 to stack2
stack2.append(stack3.pop())
moves.append("Move I from 3 to 2")

# Move F from stack2 to stack1
stack1.append(stack2.pop())
moves.append("Move F from 2 to 1")

# Move A from stack2 to stack1
stack1.append(stack2.pop())
moves.append("Move A from 2 to 1")

# Move D from stack1 to stack3
stack3.append(stack1.pop(0))
moves.append("Move D from 1 to 3")

# Print the moves
print("\n".join(moves))