# Define the initial positions of the knights
white_knights = {'B3': 'w', 'C1': 'w'}
black_knights = {'C3': 'B', 'A2': 'B'}
empty_squares = {'B1', 'C2'}

# Define a function to check if a move is valid
def is_valid_move(start, end, occupied):
    # Calculate the move difference
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    # Check if the move is an L-shape and the end square is empty
    return (dx, dy) in [(2, 1), (1, 2)] and end not in occupied

# Define a function to perform the moves
def perform_moves():
    moves = []
    occupied = set(white_knights.keys()).union(black_knights.keys()).union(empty_squares)
    
    # Move sequence to swap the knights
    # 1. w moves B3 to C2
    if is_valid_move('B3', 'C2', occupied):
        moves.append("w,B3,C2")
        occupied.remove('B3')
        occupied.add('C2')
        white_knights['C2'] = white_knights.pop('B3')
    
    # 2. B moves C3 to B3
    if is_valid_move('C3', 'B3', occupied):
        moves.append("B,C3,B3")
        occupied.remove('C3')
        occupied.add('B3')
        black_knights['B3'] = black_knights.pop('C3')
    
    # 3. w moves C1 to A2
    if is_valid_move('C1', 'A2', occupied):
        moves.append("w,C1,A2")
        occupied.remove('C1')
        occupied.add('A2')
        white_knights['A2'] = white_knights.pop('C1')
    
    # 4. B moves A2 to C1
    if is_valid_move('A2', 'C1', occupied):
        moves.append("B,A2,C1")
        occupied.remove('A2')
        occupied.add('C1')
        black_knights['C1'] = black_knights.pop('A2')
    
    return moves

# Get the sequence of moves
move_sequence = perform_moves()
print(move_sequence)