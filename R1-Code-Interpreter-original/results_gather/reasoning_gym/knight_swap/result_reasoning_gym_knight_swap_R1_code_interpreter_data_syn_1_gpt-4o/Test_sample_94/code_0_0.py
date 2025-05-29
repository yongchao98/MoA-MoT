def is_valid_move(start, end):
    # Calculate the difference in positions
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    # Check if the move is a valid knight move
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def swap_knights():
    # Initial positions
    white_knights = {'C3', 'B1'}
    black_knights = {'A3', 'C1'}
    empty_squares = {'A1', 'B3', 'D3', 'D2', 'D1'}

    # Target positions
    target_white = {'A3', 'C1'}
    target_black = {'C3', 'B1'}

    # Moves list
    moves = []

    # Simulate moves
    # White moves first
    # Try to move C3 to A3
    if is_valid_move('C3', 'A3') and 'A3' in empty_squares:
        moves.append("w,C3,A3")
        white_knights.remove('C3')
        white_knights.add('A3')
        empty_squares.remove('A3')
        empty_squares.add('C3')

    # Black moves
    # Try to move A3 to C3
    if is_valid_move('A3', 'C3') and 'C3' in empty_squares:
        moves.append("B,A3,C3")
        black_knights.remove('A3')
        black_knights.add('C3')
        empty_squares.remove('C3')
        empty_squares.add('A3')

    # Continue the process for other knights
    # Check if all knights are in target positions
    if white_knights == target_white and black_knights == target_black:
        return moves
    else:
        return "No"

result = swap_knights()
print(result)