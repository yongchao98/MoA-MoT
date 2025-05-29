def is_valid_move(start, end):
    # Calculate the difference in positions
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    # Check if the move is a valid knight move
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def find_moves():
    # Initial positions
    white_knights = {'C3', 'A2'}
    black_knights = {'B1', 'D1'}
    empty_squares = {'A1', 'B2', 'C1', 'D2', 'D3'}

    # Target positions
    target_white = {'B1', 'D1'}
    target_black = {'C3', 'A2'}

    # Moves list
    moves = []

    # Simulate moves
    # This is a simplified simulation, not a full search
    # We will try to move white knights first, then black knights
    # This is a heuristic approach and may not find a solution if it exists

    # Move white knight from C3 to B1
    if is_valid_move('C3', 'B1') and 'B1' in empty_squares:
        moves.append("w,C3,B1")
        white_knights.remove('C3')
        white_knights.add('B1')
        empty_squares.remove('B1')
        empty_squares.add('C3')

    # Move black knight from B1 to C3
    if is_valid_move('B1', 'C3') and 'C3' in empty_squares:
        moves.append("B,B1,C3")
        black_knights.remove('B1')
        black_knights.add('C3')
        empty_squares.remove('C3')
        empty_squares.add('B1')

    # Move white knight from A2 to D1
    if is_valid_move('A2', 'D1') and 'D1' in empty_squares:
        moves.append("w,A2,D1")
        white_knights.remove('A2')
        white_knights.add('D1')
        empty_squares.remove('D1')
        empty_squares.add('A2')

    # Move black knight from D1 to A2
    if is_valid_move('D1', 'A2') and 'A2' in empty_squares:
        moves.append("B,D1,A2")
        black_knights.remove('D1')
        black_knights.add('A2')
        empty_squares.remove('A2')
        empty_squares.add('D1')

    # Check if all knights are in target positions
    if white_knights == target_white and black_knights == target_black:
        return moves
    else:
        return "No"

# Execute the function and print the result
print(find_moves())