# Initial board setup
board = [
    ['.', '.', 'B', 'B', 'B', 'I'],
    ['x', '.', 'C', 'C', 'C', 'I'],
    ['A', 'A', '.', 'G', 'H', 'I'],
    ['.', '.', '.', 'G', 'H', '.'],
    ['.', '.', 'D', 'D', 'x', '.'],
    ['.', '.', 'E', 'E', 'F', 'F']
]

# Moves to solve the puzzle
moves = [
    ('I', 1),  # Move 'I' down
    ('B', -1), # Move 'B' left
    ('C', -1), # Move 'C' up
    ('G', 1),  # Move 'G' down
    ('H', 1),  # Move 'H' down
    ('A', 3)   # Move 'A' right to exit
]

# Function to apply moves
def apply_moves(board, moves):
    for car, move in moves:
        if car in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']:
            if move > 0:
                print(f"{car}+{move}", end=' ')
            else:
                print(f"{car}{move}", end=' ')
    print()

apply_moves(board, moves)