# Define the initial and goal states
initial_state = {
    'stack1': ['D', 'E', 'H', 'G'],
    'stack2': ['A', 'F', 'I'],
    'stack3': ['B', 'C']
}

goal_state = {
    'stack1': ['G', 'I'],
    'stack2': ['A', 'C', 'D', 'F', 'H'],
    'stack3': ['B', 'E']
}

# Define the moves to reach the goal state
moves = [
    ('G', 'stack1', 'stack3'),
    ('H', 'stack1', 'stack2'),
    ('E', 'stack1', 'stack3'),
    ('D', 'stack1', 'stack2'),
    ('I', 'stack2', 'stack1'),
    ('C', 'stack3', 'stack2'),
    ('G', 'stack3', 'stack1')
]

# Print the moves
for move in moves:
    print(f"Move {move[0]} from {move[1]} to {move[2]}")