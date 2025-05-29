# Initial and goal states
initial_state = {
    'stack1': ['H', 'E', 'D', 'F', 'A'],
    'stack2': ['G', 'B'],
    'stack3': ['C']
}

goal_state = {
    'stack1': ['G'],
    'stack2': ['A', 'B', 'D', 'E'],
    'stack3': ['C', 'F', 'H']
}

# Function to print the moves
def print_moves(moves):
    for move in moves:
        print(f"Move {move[0]} from {move[1]} to {move[2]}")

# Plan of moves
moves = [
    ('A', 'stack1', 'stack2'),
    ('F', 'stack1', 'stack3'),
    ('D', 'stack1', 'stack2'),
    ('E', 'stack1', 'stack2'),
    ('H', 'stack1', 'stack3'),
    ('G', 'stack2', 'stack1'),
    ('B', 'stack2', 'stack2'),
    ('A', 'stack2', 'stack2'),
    ('F', 'stack3', 'stack3'),
    ('H', 'stack3', 'stack3')
]

# Execute the moves
print_moves(moves)