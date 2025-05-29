def is_valid_move(start, end):
    # Calculate the difference in positions
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    # Check if the move is a valid knight move
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def find_moves():
    # Initial positions
    white_knights = {'B3', 'C3'}
    black_knights = {'D2', 'C1'}
    empty_squares = {'A3', 'A1', 'B1', 'D1'}

    # Target positions
    target_white = {'D2', 'C1'}
    target_black = {'B3', 'C3'}

    # Moves list
    moves = []

    # Simulate moves
    # Move white knight from B3 to D1
    if is_valid_move('B3', 'D1'):
        moves.append("w,B3,D1")
        white_knights.remove('B3')
        white_knights.add('D1')
        empty_squares.remove('D1')
        empty_squares.add('B3')

    # Move black knight from C1 to B3
    if is_valid_move('C1', 'B3'):
        moves.append("B,C1,B3")
        black_knights.remove('C1')
        black_knights.add('B3')
        empty_squares.remove('B3')
        empty_squares.add('C1')

    # Move white knight from C3 to C1
    if is_valid_move('C3', 'C1'):
        moves.append("w,C3,C1")
        white_knights.remove('C3')
        white_knights.add('C1')
        empty_squares.remove('C1')
        empty_squares.add('C3')

    # Move black knight from D2 to C3
    if is_valid_move('D2', 'C3'):
        moves.append("B,D2,C3")
        black_knights.remove('D2')
        black_knights.add('C3')
        empty_squares.remove('C3')
        empty_squares.add('D2')

    # Move white knight from D1 to D2
    if is_valid_move('D1', 'D2'):
        moves.append("w,D1,D2")
        white_knights.remove('D1')
        white_knights.add('D2')
        empty_squares.remove('D2')
        empty_squares.add('D1')

    # Check if all knights are in target positions
    if white_knights == target_white and black_knights == target_black:
        return moves
    else:
        return "No"

# Execute the function and print the result
result = find_moves()
print(result)