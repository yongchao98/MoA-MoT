def is_valid_move(start, end, occupied):
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    return (dx, dy) in [(2, 1), (1, 2)] and end not in occupied

def swap_knights():
    # Initial positions
    white_knights = {'B1', 'A2'}
    black_knights = {'C3', 'B2'}
    empty_squares = {'B3', 'C1'}
    
    # Moves list
    moves = []
    
    # White Move 1: B1 to C3
    if is_valid_move('B1', 'C3', white_knights | black_knights):
        moves.append("w,B1,C3")
        white_knights.remove('B1')
        white_knights.add('C3')
        empty_squares.remove('C3')
        empty_squares.add('B1')
    
    # Black Move 1: C3 to B1
    if is_valid_move('C3', 'B1', white_knights | black_knights):
        moves.append("B,C3,B1")
        black_knights.remove('C3')
        black_knights.add('B1')
        empty_squares.remove('B1')
        empty_squares.add('C3')
    
    # White Move 2: A2 to B2
    if is_valid_move('A2', 'B2', white_knights | black_knights):
        moves.append("w,A2,B2")
        white_knights.remove('A2')
        white_knights.add('B2')
        empty_squares.remove('B2')
        empty_squares.add('A2')
    
    # Black Move 2: B2 to A2
    if is_valid_move('B2', 'A2', white_knights | black_knights):
        moves.append("B,B2,A2")
        black_knights.remove('B2')
        black_knights.add('A2')
        empty_squares.remove('A2')
        empty_squares.add('B2')
    
    # Check if the swap is complete
    if white_knights == {'C3', 'B2'} and black_knights == {'B1', 'A2'}:
        return moves
    else:
        return "No"

# Execute the function and print the result
print(swap_knights())