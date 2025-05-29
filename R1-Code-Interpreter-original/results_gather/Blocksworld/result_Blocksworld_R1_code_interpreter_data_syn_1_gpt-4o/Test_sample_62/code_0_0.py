# Initial state
stack1 = ['B']
stack2 = ['G', 'A', 'F', 'E']
stack3 = ['D', 'C', 'H']

# Moves to achieve the goal state
moves = []

# Move C from stack3 to stack2
stack2.append(stack3.pop())
moves.append("Move C from 3 to 2")

# Move H from stack3 to stack2 temporarily
stack2.append(stack3.pop())
moves.append("Move H from 3 to 2")

# Move E from stack2 to stack3
stack3.append(stack2.pop())
moves.append("Move E from 2 to 3")

# Move F from stack2 to stack3
stack3.append(stack2.pop())
moves.append("Move F from 2 to 3")

# Move A from stack2 to stack3
stack3.append(stack2.pop())
moves.append("Move A from 2 to 3")

# Move H back from stack2 to stack3
stack3.append(stack2.pop())
moves.append("Move H from 2 to 3")

# Print the moves
for move in moves:
    print(move)