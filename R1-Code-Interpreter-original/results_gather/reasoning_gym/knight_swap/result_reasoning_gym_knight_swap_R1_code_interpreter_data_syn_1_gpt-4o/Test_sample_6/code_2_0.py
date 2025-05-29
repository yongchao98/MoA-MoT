# Define the initial positions of the knights
white_knights = {'D3': 'w', 'D1': 'w'}
black_knights = {'C3': 'B', 'C1': 'B'}
empty_squares = ['A2', 'B1']

# Define the moves
moves = [
    ("w", "D3", "B2"),  # White knight moves from D3 to B2
    ("B", "C3", "A2"),  # Black knight moves from C3 to A2
    ("w", "D1", "C3"),  # White knight moves from D1 to C3
    ("B", "C1", "D3"),  # Black knight moves from C1 to D3
    ("w", "B2", "C1"),  # White knight moves from B2 to C1
    ("B", "A2", "D1")   # Black knight moves from A2 to D1
]

# Function to apply moves and check the final positions
def apply_moves(moves, white_knights, black_knights):
    for move in moves:
        color, from_pos, to_pos = move
        if color == 'w':
            # Move white knight
            white_knights[to_pos] = white_knights.pop(from_pos)
        else:
            # Move black knight
            black_knights[to_pos] = black_knights.pop(from_pos)
    
    # Check if the swap is successful
    final_white_positions = set(black_knights.keys())
    final_black_positions = set(white_knights.keys())
    
    return final_white_positions == set(['C3', 'C1']) and final_black_positions == set(['D3', 'D1'])

# Check if the moves achieve the swap
is_possible = apply_moves(moves, white_knights, black_knights)
print("Yes" if is_possible else "No")
if is_possible:
    print(moves)