def is_valid_move(start, end):
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def find_moves():
    # Initial positions
    white_knights = {'A1': 'w', 'C1': 'w'}
    black_knights = {'B2': 'B', 'D2': 'B'}
    empty_squares = {'A3', 'B3', 'C3', 'D3'}

    # Possible moves
    moves = []

    # Black moves first
    # Move B from B2 to A3
    if is_valid_move('B2', 'A3'):
        moves.append('B,B2,A3')
        black_knights.pop('B2')
        black_knights['A3'] = 'B'
        empty_squares.remove('A3')
        empty_squares.add('B2')

    # White moves
    # Move w from C1 to B3
    if is_valid_move('C1', 'B3'):
        moves.append('w,C1,B3')
        white_knights.pop('C1')
        white_knights['B3'] = 'w'
        empty_squares.remove('B3')
        empty_squares.add('C1')

    # Black moves
    # Move B from D2 to C3
    if is_valid_move('D2', 'C3'):
        moves.append('B,D2,C3')
        black_knights.pop('D2')
        black_knights['C3'] = 'B'
        empty_squares.remove('C3')
        empty_squares.add('D2')

    # White moves
    # Move w from A1 to D2
    if is_valid_move('A1', 'D2'):
        moves.append('w,A1,D2')
        white_knights.pop('A1')
        white_knights['D2'] = 'w'
        empty_squares.remove('D2')
        empty_squares.add('A1')

    # Black moves
    # Move B from A3 to C1
    if is_valid_move('A3', 'C1'):
        moves.append('B,A3,C1')
        black_knights.pop('A3')
        black_knights['C1'] = 'B'
        empty_squares.remove('C1')
        empty_squares.add('A3')

    # White moves
    # Move w from B3 to B2
    if is_valid_move('B3', 'B2'):
        moves.append('w,B3,B2')
        white_knights.pop('B3')
        white_knights['B2'] = 'w'
        empty_squares.remove('B2')
        empty_squares.add('B3')

    # Check if all knights are in their target positions
    if white_knights == {'B2': 'w', 'D2': 'w'} and black_knights == {'A1': 'B', 'C1': 'B'}:
        return moves
    else:
        return "No"

# Execute the function and print the result
print(find_moves())