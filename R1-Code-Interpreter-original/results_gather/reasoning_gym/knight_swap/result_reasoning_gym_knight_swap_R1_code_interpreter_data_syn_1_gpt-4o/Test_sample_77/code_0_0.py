def is_valid_move(start, end):
    # Calculate the difference in coordinates
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    # Check if the move is a valid knight move
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def find_moves():
    # Initial positions
    white_knights = ['B1', 'D2']
    black_knights = ['A3', 'C2']
    empty_squares = ['A1', 'C1', 'B3']

    # Target positions
    target_white = ['A3', 'C2']
    target_black = ['B1', 'D2']

    # Moves list
    moves = []

    # Simulate moves
    # Move 1: w B1 to C3
    if is_valid_move('B1', 'C3') and 'C3' in empty_squares:
        moves.append('w,B1,C3')
        white_knights.remove('B1')
        white_knights.append('C3')
        empty_squares.remove('C3')
        empty_squares.append('B1')

    # Move 2: B A3 to B1
    if is_valid_move('A3', 'B1') and 'B1' in empty_squares:
        moves.append('B,A3,B1')
        black_knights.remove('A3')
        black_knights.append('B1')
        empty_squares.remove('B1')
        empty_squares.append('A3')

    # Move 3: w D2 to A3
    if is_valid_move('D2', 'A3') and 'A3' in empty_squares:
        moves.append('w,D2,A3')
        white_knights.remove('D2')
        white_knights.append('A3')
        empty_squares.remove('A3')
        empty_squares.append('D2')

    # Move 4: B C2 to D2
    if is_valid_move('C2', 'D2') and 'D2' in empty_squares:
        moves.append('B,C2,D2')
        black_knights.remove('C2')
        black_knights.append('D2')
        empty_squares.remove('D2')
        empty_squares.append('C2')

    # Move 5: w C3 to C2
    if is_valid_move('C3', 'C2') and 'C2' in empty_squares:
        moves.append('w,C3,C2')
        white_knights.remove('C3')
        white_knights.append('C2')
        empty_squares.remove('C2')
        empty_squares.append('C3')

    # Check if all knights are in target positions
    if set(white_knights) == set(target_white) and set(black_knights) == set(target_black):
        return moves
    else:
        return "No"

# Execute the function and print the result
print(find_moves())