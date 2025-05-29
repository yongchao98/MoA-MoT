# Initial board setup
board = [
    ['.', '.', '.', 'E', '.', '.'],
    ['.', '.', 'C', 'E', '.', '.'],
    ['A', 'A', 'C', 'E', '.', '.'],
    ['.', '.', '.', '.', '.', '.'],
    ['.', '.', 'D', '.', '.', '.'],
    ['.', '.', 'D', 'B', 'B', 'x']
]

# Define the moves
moves = []

# Move C up by 1
moves.append('C-1')

# Move D down by 1
moves.append('D+1')

# Move C right by 2
moves.append('C+2')

# Move A right by 2 to reach the exit
moves.append('A+2')

# Print the sequence of moves
print(' '.join(moves))