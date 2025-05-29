def is_valid_move(start, end):
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def find_moves():
    # Initial positions
    white_knights = {'A3': 'w', 'D2': 'w'}
    black_knights = {'B3': 'B', 'C2': 'B'}
    empty_squares = {'A1', 'B1', 'C1', 'C3', 'B4', 'D4', 'E1', 'E4'}

    # Target positions
    target_white = {'B3', 'C2'}
    target_black = {'A3', 'D2'}

    # Moves list
    moves = []

    # Move white knight from A3 to C2
    if is_valid_move('A3', 'C2') and 'C2' in empty_squares:
        moves.append("w,A3,C2")
        white_knights.pop('A3')
        white_knights['C2'] = 'w'
        empty_squares.remove('C2')
        empty_squares.add('A3')

    # Move black knight from C2 to A3
    if is_valid_move('C2', 'A3') and 'A3' in empty_squares:
        moves.append("B,C2,A3")
        black_knights.pop('C2')
        black_knights['A3'] = 'B'
        empty_squares.remove('A3')
        empty_squares.add('C2')

    # Move white knight from D2 to B3
    if is_valid_move('D2', 'B3') and 'B3' in empty_squares:
        moves.append("w,D2,B3")
        white_knights.pop('D2')
        white_knights['B3'] = 'w'
        empty_squares.remove('B3')
        empty_squares.add('D2')

    # Move black knight from B3 to D2
    if is_valid_move('B3', 'D2') and 'D2' in empty_squares:
        moves.append("B,B3,D2")
        black_knights.pop('B3')
        black_knights['D2'] = 'B'
        empty_squares.remove('D2')
        empty_squares.add('B3')

    # Check if the swap is complete
    if set(white_knights.keys()) == target_white and set(black_knights.keys()) == target_black:
        return moves
    else:
        return "No"

# Execute the function and print the result
print(find_moves())