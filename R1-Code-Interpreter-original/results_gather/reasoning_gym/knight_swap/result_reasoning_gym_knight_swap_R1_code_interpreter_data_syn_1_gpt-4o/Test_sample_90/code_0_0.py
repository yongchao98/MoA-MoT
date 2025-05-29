def is_valid_move(start, end, board):
    # Calculate the move difference
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    # Check if it's a valid knight move
    return (dx == 2 and dy == 1) or (dx == 1 and dy == 2)

def swap_knights():
    # Initial positions
    white_knights = ['C1', 'D1']
    black_knights = ['A3', 'A2']
    empty_squares = ['B1', 'B2', 'B3', 'B4', 'C2', 'C3', 'C4', 'D2', 'D3']
    
    # Target positions
    target_white = ['A3', 'A2']
    target_black = ['C1', 'D1']
    
    # Moves list
    moves = []
    
    # Simulate moves
    # White moves first
    # Try moving C1 to A2
    if is_valid_move('C1', 'A2', empty_squares):
        moves.append("w,C1,A2")
        white_knights.remove('C1')
        white_knights.append('A2')
        empty_squares.remove('A2')
        empty_squares.append('C1')
    
    # Black moves
    # Try moving A3 to C1
    if is_valid_move('A3', 'C1', empty_squares):
        moves.append("B,A3,C1")
        black_knights.remove('A3')
        black_knights.append('C1')
        empty_squares.remove('C1')
        empty_squares.append('A3')
    
    # White moves
    # Try moving D1 to A3
    if is_valid_move('D1', 'A3', empty_squares):
        moves.append("w,D1,A3")
        white_knights.remove('D1')
        white_knights.append('A3')
        empty_squares.remove('A3')
        empty_squares.append('D1')
    
    # Black moves
    # Try moving A2 to D1
    if is_valid_move('A2', 'D1', empty_squares):
        moves.append("B,A2,D1")
        black_knights.remove('A2')
        black_knights.append('D1')
        empty_squares.remove('D1')
        empty_squares.append('A2')
    
    # Check if the swap is complete
    if set(white_knights) == set(target_white) and set(black_knights) == set(target_black):
        return moves
    else:
        return "No"

# Execute the function and print the result
print(swap_knights())