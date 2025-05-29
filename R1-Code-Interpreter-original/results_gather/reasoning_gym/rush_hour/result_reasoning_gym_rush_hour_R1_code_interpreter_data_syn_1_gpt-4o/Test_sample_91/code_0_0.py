# Initial board setup
board = [
    "BBCCC.",
    "....HI",
    "..AAHI",
    "...GDD",
    "...G..",
    ".xEEFF"
]

# Function to find the position of a car
def find_car_position(board, car):
    positions = []
    for r, row in enumerate(board):
        for c, cell in enumerate(row):
            if cell == car:
                positions.append((r, c))
    return positions

# Function to generate moves
def generate_moves():
    moves = []
    # Move HI down
    moves.append("H+1 I+1")
    # Move G down
    moves.append("G+1")
    # Move DD down
    moves.append("D+1 D+1")
    # Move AA to the right
    moves.append("A+2")
    return moves

# Generate the sequence of moves
moves = generate_moves()

# Print the moves in the required format
print("<<<" + " ".join(moves) + ">>>")