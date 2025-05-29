def is_valid_move(start, end):
    # Calculate the difference in positions
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    # Check if the move is a valid knight move
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def swap_knights():
    # Initial positions
    white_knights = ['C3', 'B1']
    black_knights = ['A3', 'D1']
    empty_squares = ['B3', 'A1', 'D3', 'C1']

    # Target positions
    target_white = ['A3', 'D1']
    target_black = ['C3', 'B1']

    # Moves list
    moves = []

    # Simulate moves
    # Move 1: White knight from C3 to D1
    if is_valid_move('C3', 'D1') and 'D1' in empty_squares:
        white_knights[0] = 'D1'
        empty_squares.remove('D1')
        empty_squares.append('C3')
        moves.append("w,C3,D1")

    # Move 2: Black knight from A3 to C3
    if is_valid_move('A3', 'C3') and 'C3' in empty_squares:
        black_knights[0] = 'C3'
        empty_squares.remove('C3')
        empty_squares.append('A3')
        moves.append("B,A3,C3")

    # Move 3: White knight from B1 to A3
    if is_valid_move('B1', 'A3') and 'A3' in empty_squares:
        white_knights[1] = 'A3'
        empty_squares.remove('A3')
        empty_squares.append('B1')
        moves.append("w,B1,A3")

    # Move 4: Black knight from D1 to B1
    if is_valid_move('D1', 'B1') and 'B1' in empty_squares:
        black_knights[1] = 'B1'
        empty_squares.remove('B1')
        empty_squares.append('D1')
        moves.append("B,D1,B1")

    # Check if the knights are in their target positions
    if sorted(white_knights) == sorted(target_white) and sorted(black_knights) == sorted(target_black):
        return moves
    else:
        return "No"

# Execute the function and print the result
print(swap_knights())