# Define the moves to transform the initial state to the goal state
moves = [
    "Move D from 4 to 2",
    "Move H from 4 to 1",
    "Move C from 4 to 3",
    "Move B from 2 to 1",
    "Move F from 2 to 1",
    "Move E from 3 to 1",
    "Move I from 1 to 2",
    "Move G from 1 to 3",
    "Move J from 3 to 3",  # J is already in the correct stack
    "Move K from 1 to 3"
]

# Print the moves
for move in moves:
    print(move)