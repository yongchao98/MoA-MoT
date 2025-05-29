# Initial state of the stacks
stack1 = ['B', 'K', 'J', 'C']
stack2 = ['F', 'H']
stack3 = ['D', 'E', 'G', 'I', 'A']

# Goal state of the stacks
goal_stack1 = ['C', 'D', 'F', 'I', 'K']
goal_stack2 = ['B', 'H', 'J']
goal_stack3 = ['A', 'E', 'G']

# Function to print the move
def move(block, source, destination):
    print(f"Move {block} from {source} to {destination}")

# Simulate the moves
# Move A from stack3 to stack2
stack2.append(stack3.pop())
move('A', 'stack3', 'stack2')

# Move I from stack3 to stack1
stack1.append(stack3.pop())
move('I', 'stack3', 'stack1')

# Move G from stack3 to stack2
stack2.append(stack3.pop())
move('G', 'stack3', 'stack2')

# Move E from stack3 to stack1
stack1.append(stack3.pop())
move('E', 'stack3', 'stack1')

# Move D from stack3 to stack1
stack1.append(stack3.pop())
move('D', 'stack3', 'stack1')

# Move H from stack2 to stack3
stack3.append(stack2.pop())
move('H', 'stack2', 'stack3')

# Move F from stack2 to stack1
stack1.append(stack2.pop())
move('F', 'stack2', 'stack1')

# Move J from stack1 to stack2
stack2.append(stack1.pop())
move('J', 'stack1', 'stack2')

# Move K from stack1 to stack3
stack3.append(stack1.pop())
move('K', 'stack1', 'stack3')

# Move B from stack1 to stack2
stack2.append(stack1.pop())
move('B', 'stack1', 'stack2')

# Move K from stack3 to stack1
stack1.append(stack3.pop())
move('K', 'stack3', 'stack1')

# Move H from stack3 to stack2
stack2.append(stack3.pop())
move('H', 'stack3', 'stack2')

# Move G from stack2 to stack3
stack3.append(stack2.pop())
move('G', 'stack2', 'stack3')

# Move A from stack2 to stack3
stack3.append(stack2.pop())
move('A', 'stack2', 'stack3')