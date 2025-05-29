# Initial and goal states
initial_state = {
    'stack1': ['B', 'G'],
    'stack2': ['F', 'E', 'A', 'C'],
    'stack3': ['H', 'D']
}

goal_state = {
    'stack1': ['A', 'B', 'F'],
    'stack2': ['D', 'E'],
    'stack3': ['C', 'G', 'H']
}

# Function to print the moves
def print_moves(moves):
    for move in moves:
        print(move)

# Plan of moves to reach the goal state
moves = [
    "Move G from stack1 to stack3",
    "Move B from stack1 to stack2",
    "Move C from stack2 to stack3",
    "Move A from stack2 to stack1",
    "Move F from stack2 to stack1",
    "Move D from stack3 to stack2",
    "Move B from stack2 to stack1"
]

# Print the moves
print_moves(moves)