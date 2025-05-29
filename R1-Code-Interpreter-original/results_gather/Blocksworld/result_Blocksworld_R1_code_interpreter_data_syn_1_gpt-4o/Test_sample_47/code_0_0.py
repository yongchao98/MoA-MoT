# Initial and goal states
initial_state = {
    'stack1': ['D', 'G', 'E', 'C'],
    'stack2': ['A'],
    'stack3': ['F'],
    'stack4': ['B']
}

goal_state = {
    'stack1': ['D', 'G'],
    'stack2': ['A', 'C'],
    'stack3': ['B', 'E', 'F']
}

# Moves to achieve the goal state
moves = []

# Move C from stack1 to stack2
moves.append("Move C from stack1 to stack2")
initial_state['stack2'].append(initial_state['stack1'].pop())

# Move E from stack1 to stack3
moves.append("Move E from stack1 to stack3")
initial_state['stack3'].append(initial_state['stack1'].pop())

# Move B from stack4 to stack3
moves.append("Move B from stack4 to stack3")
initial_state['stack3'].append(initial_state['stack4'].pop())

# Print the moves
print("\n".join(moves))