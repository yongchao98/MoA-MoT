def is_valid_move(start, end):
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def swap_knights():
    # Initial positions
    white_knights = {'A2', 'B1'}
    black_knights = {'C3', 'C2'}
    empty_squares = {'A1', 'A3', 'B3', 'B4', 'C1', 'D1', 'D2', 'D3', 'D4'}
    
    # Target positions
    target_white = {'C3', 'C2'}
    target_black = {'A2', 'B1'}
    
    # Moves list
    moves = []
    
    # Simulate moves
    # Black moves first
    # Try to move black knights to target positions
    if is_valid_move('C3', 'B1') and 'B1' in empty_squares:
        moves.append("B,C3,B1")
        black_knights.remove('C3')
        black_knights.add('B1')
        empty_squares.remove('B1')
        empty_squares.add('C3')
    
    if is_valid_move('C2', 'A2') and 'A2' in empty_squares:
        moves.append("B,C2,A2")
        black_knights.remove('C2')
        black_knights.add('A2')
        empty_squares.remove('A2')
        empty_squares.add('C2')
    
    # Now move white knights to target positions
    if is_valid_move('A2', 'C3') and 'C3' in empty_squares:
        moves.append("w,A2,C3")
        white_knights.remove('A2')
        white_knights.add('C3')
        empty_squares.remove('C3')
        empty_squares.add('A2')
    
    if is_valid_move('B1', 'C2') and 'C2' in empty_squares:
        moves.append("w,B1,C2")
        white_knights.remove('B1')
        white_knights.add('C2')
        empty_squares.remove('C2')
        empty_squares.add('B1')
    
    # Check if all knights are in target positions
    if white_knights == target_white and black_knights == target_black:
        return moves
    else:
        return "No"

# Execute the function and print the result
print(swap_knights())