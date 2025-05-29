def is_valid_move(start, end):
    # Calculate the difference in positions
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    # Check if the move is a valid knight move
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def find_moves():
    # Initial positions
    white_knights = ['A3', 'D1']
    black_knights = ['B1', 'C1']
    empty_squares = ['B3', 'A1']

    # Target positions
    target_white = ['B1', 'C1']
    target_black = ['A3', 'D1']

    # Moves list
    moves = []

    # Simulate moves
    # Move white knight from A3 to B1
    if is_valid_move('A3', 'B1') and 'B1' in empty_squares:
        moves.append("w,A3,B1")
        white_knights.remove('A3')
        white_knights.append('B1')
        empty_squares.remove('B1')
        empty_squares.append('A3')

    # Move black knight from B1 to A3
    if is_valid_move('B1', 'A3') and 'A3' in empty_squares:
        moves.append("B,B1,A3")
        black_knights.remove('B1')
        black_knights.append('A3')
        empty_squares.remove('A3')
        empty_squares.append('B1')

    # Move white knight from D1 to C1
    if is_valid_move('D1', 'C1') and 'C1' in empty_squares:
        moves.append("w,D1,C1")
        white_knights.remove('D1')
        white_knights.append('C1')
        empty_squares.remove('C1')
        empty_squares.append('D1')

    # Move black knight from C1 to D1
    if is_valid_move('C1', 'D1') and 'D1' in empty_squares:
        moves.append("B,C1,D1")
        black_knights.remove('C1')
        black_knights.append('D1')
        empty_squares.remove('D1')
        empty_squares.append('C1')

    # Check if the target configuration is reached
    if set(white_knights) == set(target_white) and set(black_knights) == set(target_black):
        return moves
    else:
        return "No"

# Execute the function and print the result
result = find_moves()
print(result)