# Define the initial board state
board = [
    "BBCC.x",
    ".DDEEL",
    ".AA..L",
    "xFFK..",
    "I.JKGG",
    "I.JHH."
]

# Define the moves
moves = [
    "L+1",  # Move L down
    "K+1",  # Move K right
    "AA+3"  # Move AA to the exit
]

# Function to apply moves and print the result
def apply_moves(board, moves):
    for move in moves:
        car, direction = move[:-2], move[-2:]
        print(f"Move {car} {direction}")

apply_moves(board, moves)