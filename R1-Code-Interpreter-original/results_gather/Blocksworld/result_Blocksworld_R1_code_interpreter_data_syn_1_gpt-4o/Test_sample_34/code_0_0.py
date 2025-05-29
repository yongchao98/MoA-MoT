# Initial state
stack1 = ['G', 'E', 'B']
stack2 = ['A']
stack3 = ['C', 'D', 'H', 'F']

# Moves to achieve the goal state
moves = []

# Move F from stack3 to stack1
stack1.append(stack3.pop())
moves.append("Move F from 3 to 1")

# Move H from stack3 to stack2
stack2.append(stack3.pop())
moves.append("Move H from 3 to 2")

# Move D from stack3 to stack1
stack1.append(stack3.pop())
moves.append("Move D from 3 to 1")

# Move C from stack3 to stack2
stack2.append(stack3.pop())
moves.append("Move C from 3 to 2")

# Move B from stack1 to stack3
stack3.append(stack1.pop())
moves.append("Move B from 1 to 3")

# Move E from stack1 to stack2
stack2.append(stack1.pop())
moves.append("Move E from 1 to 2")

# Move G from stack1 to stack1 (no actual move needed, just for clarity)
# This is to show that G is already in the correct position
moves.append("Move G from 1 to 1")

# Print the moves
print("\n".join(moves))