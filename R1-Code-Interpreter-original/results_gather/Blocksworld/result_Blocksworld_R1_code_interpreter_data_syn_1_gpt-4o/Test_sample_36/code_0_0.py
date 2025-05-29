# Initial state
stack1 = ['D', 'E']
stack2 = ['G', 'I', 'F']
stack3 = ['C', 'B', 'H', 'A']

# Moves to reach the goal state
moves = []

# Move H from stack3 to stack2
stack2.append(stack3.pop())
moves.append("Move H from 3 to 2")

# Move B from stack3 to stack1
stack1.append(stack3.pop())
moves.append("Move B from 3 to 1")

# Move C from stack3 to stack1
stack1.append(stack3.pop())
moves.append("Move C from 3 to 1")

# Move D from stack1 to stack3
stack3.append(stack1.pop(0))
moves.append("Move D from 1 to 3")

# Move F from stack2 to stack1
stack1.append(stack2.pop())
moves.append("Move F from 2 to 1")

# Move I from stack2 to stack1
stack1.append(stack2.pop())
moves.append("Move I from 2 to 1")

# Print the moves
print("\n".join(moves))