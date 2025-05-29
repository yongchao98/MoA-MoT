def is_valid_move(start, end, board):
    # Calculate the move difference
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    # Check if it's a valid knight move
    return (dx, dy) in [(2, 1), (1, 2)] and board[end] == '.'

def swap_knights():
    # Initial positions
    board = {
        'A3': '.', 'B3': '.', 'C3': ' ', 'D3': '.',
        'A2': ' ', 'B2': 'w', 'C2': 'w', 'D2': '.',
        'A1': ' ', 'B1': '.', 'C1': 'B', 'D1': 'B'
    }
    
    # Target positions
    target = {
        'B2': 'B', 'C2': 'B',
        'C1': 'w', 'D1': 'w'
    }
    
    # Moves list
    moves = []
    
    # Simulate moves
    def move_knight(color, start, end):
        if is_valid_move(start, end, board):
            board[end] = color
            board[start] = '.'
            moves.append(f"{color},{start},{end}")
            return True
        return False
    
    # Attempt to swap positions
    # This is a simplified attempt to find a solution
    # White moves
    if move_knight('w', 'B2', 'D1'):
        # Black moves
        if move_knight('B', 'C1', 'B2'):
            # White moves
            if move_knight('w', 'C2', 'C1'):
                # Black moves
                if move_knight('B', 'D1', 'C2'):
                    # Check if target is reached
                    if all(board[pos] == target[pos] for pos in target):
                        return moves
    return "No"

print(swap_knights())