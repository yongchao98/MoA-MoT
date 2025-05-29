def is_valid_move(start, end, occupied):
    # Calculate the difference in positions
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    # Check if the move is a valid knight move and the end position is not occupied
    return (dx, dy) in [(2, 1), (1, 2)] and end not in occupied

def swap_knights():
    # Initial positions
    white_knights = ['B1', 'D1']
    black_knights = ['C3', 'B2']
    empty_squares = ['A3', 'A2', 'C1']
    
    # Target positions
    target_white = ['C3', 'B2']
    target_black = ['B1', 'D1']
    
    # Moves list
    moves = []
    
    # Simulate moves
    # First move: w from B1 to C3
    if is_valid_move('B1', 'C3', black_knights + empty_squares):
        white_knights.remove('B1')
        white_knights.append('C3')
        empty_squares.remove('C3')
        empty_squares.append('B1')
        moves.append("w,B1,C3")
    
    # Second move: B from C3 to B1
    if is_valid_move('C3', 'B1', white_knights + empty_squares):
        black_knights.remove('C3')
        black_knights.append('B1')
        empty_squares.remove('B1')
        empty_squares.append('C3')
        moves.append("B,C3,B1")
    
    # Third move: w from D1 to B2
    if is_valid_move('D1', 'B2', black_knights + empty_squares):
        white_knights.remove('D1')
        white_knights.append('B2')
        empty_squares.remove('B2')
        empty_squares.append('D1')
        moves.append("w,D1,B2")
    
    # Fourth move: B from B2 to D1
    if is_valid_move('B2', 'D1', white_knights + empty_squares):
        black_knights.remove('B2')
        black_knights.append('D1')
        empty_squares.remove('D1')
        empty_squares.append('B2')
        moves.append("B,B2,D1")
    
    # Check if the swap is complete
    if set(white_knights) == set(target_white) and set(black_knights) == set(target_black):
        return moves
    else:
        return "No"

# Execute the function and print the result
print(swap_knights())