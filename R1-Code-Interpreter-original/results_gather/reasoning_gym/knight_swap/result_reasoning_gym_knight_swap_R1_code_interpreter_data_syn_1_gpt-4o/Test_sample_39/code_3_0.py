def is_valid_move(start, end, board):
    dx = abs(ord(end[0]) - ord(start[0]))
    dy = abs(int(end[1]) - int(start[1]))
    return (dx, dy) in [(2, 1), (1, 2)] and board[end] == '.'

def swap_knights():
    board = {
        'A1': '.', 'B1': 'B', 'C1': '.', 'D1': '.',
        'A2': 'B', 'B2': '.', 'C2': 'w', 'D2': 'w',
        'A3': '.', 'B3': 'B', 'C3': '.', 'D3': '.',
        'A4': '.', 'B4': '.', 'C4': '.', 'D4': '.'
    }
    
    target_positions = {
        'w': ['A2', 'B1', 'B3'],
        'B': ['C2', 'D2']
    }
    
    moves = []
    knight_positions = {'w': ['C2', 'D2'], 'B': ['A2', 'B1', 'B3']}
    
    def backtrack(turn):
        if all(board[pos] == 'w' for pos in target_positions['w']) and \
           all(board[pos] == 'B' for pos in target_positions['B']):
            return True
        
        color = 'B' if turn % 2 == 0 else 'w'
        for i, start in enumerate(knight_positions[color]):
            for dx, dy in [(2, 1), (1, 2), (-1, 2), (-2, 1), (-2, -1), (-1, -2), (1, -2), (2, -1)]:
                end = chr(ord(start[0]) + dx) + str(int(start[1]) + dy)
                if end in board and is_valid_move(start, end, board):
                    board[end] = color
                    board[start] = '.'
                    knight_positions[color][i] = end
                    moves.append(f"{color},{start},{end}")
                    
                    if backtrack(turn + 1):
                        return True
                    
                    # Backtrack
                    board[start] = color
                    board[end] = '.'
                    knight_positions[color][i] = start
                    moves.pop()
        
        return False
    
    if backtrack(0):
        return moves
    else:
        return "No"

print(swap_knights())