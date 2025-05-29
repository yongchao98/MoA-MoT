def swap_knights():
    # Initial positions
    board = {
        'A3': '.', 'B3': '.', 'C3': ' ', 'D3': '.',
        'A2': ' ', 'B2': 'w', 'C2': 'w', 'D2': '.',
        'A1': ' ', 'B1': '.', 'C1': 'B', 'D1': 'B'
    }
    
    # Moves list
    moves = []
    
    # Manually simulate a sequence of moves
    # White moves
    moves.append("w,B2,D1")
    board['B2'], board['D1'] = '.', 'w'
    
    # Black moves
    moves.append("B,C1,B3")
    board['C1'], board['B3'] = '.', 'B'
    
    # White moves
    moves.append("w,C2,C1")
    board['C2'], board['C1'] = '.', 'w'
    
    # Black moves
    moves.append("B,D1,C2")
    board['D1'], board['C2'] = '.', 'B'
    
    # Check if target is reached
    target = {
        'B2': 'B', 'C2': 'B',
        'C1': 'w', 'D1': 'w'
    }
    
    if all(board[pos] == target[pos] for pos in target):
        return moves
    else:
        return "No"

print(swap_knights())